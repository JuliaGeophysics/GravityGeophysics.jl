# 3D mesh design for gravity data: core, padding, layers, topography, start model and inversion mask
# Author: @pankajkmishra
# The MT mesh designer with the rules of potential-field practice: a uniform core over the stations with
# cells of about half the median station spacing (an odd count, one cell beyond the stations, so a regular
# grid of stations sits on cell centres, off the mesh planes); lateral padding growing by a factor until it
# is at least as wide as the core (potential fields decay as distance⁻², but the padding still carries the
# regional field of the edge structure); regular layers from the highest ground of the core down to the
# depth of investigation below its lowest ground (default half the core width), then layers growing by the
# vertical factor until they are as thick as the regular block. There is no skin depth. The ground comes
# from a GeoTIFF DEM or a function of (x, y) (e.g. GPS height minus radar altimeter of an airborne survey);
# station heights are kept as measured, and stations must lie on or above the ground (gz is finite on it, so
# ground stations need no clearance): the columns where one does not are lowered. Writes the start model and an
# inversion mask (1 earth, air NaN) instead of a covariance, with z = 0 at the highest ground (the top of the highest
# earth cell, the model datum, see Model.jl); the data file keeps z = −elevation

using ArchGDAL
using Printf

const _TOPO_CACHE = Dict{String, Any}()

_snap(v, r) = collect(r)[argmin(abs.(collect(r) .- v))]

#----- find the typical station spacing -----
function nearest_neighbour_spacing(sx, sy)
    n = length(sx)
    n < 2 && return (median = 100.0, min = 100.0)
    dmin = fill(Inf, n)
    @inbounds for i in 1:n, j in 1:n
        i == j && continue
        d = hypot(sx[i] - sx[j], sy[i] - sy[j])
        d < dmin[i] && (dmin[i] = d)
    end
    dmin = filter(isfinite, dmin)
    return (median = median(dmin), min = minimum(dmin))
end

#----- count how many stations fall in each core cell -----
function site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)
    counts = Dict{Tuple{Int, Int}, Int}()
    nx = length(x_edges) - 1
    ny = length(y_edges) - 1
    outside = 0
    for k in eachindex(sx)
        ix = clamp(searchsortedlast(x_edges, sx[k]), 1, nx)
        iy = clamp(searchsortedlast(y_edges, sy[k]), 1, ny)
        (ix in ix_core && iy in iy_core) || (outside += 1)
        counts[(ix, iy)] = get(counts, (ix, iy), 0) + 1
    end
    maxper = isempty(counts) ? 0 : maximum(values(counts))
    incore = count(p -> first(p)[1] in ix_core && first(p)[2] in iy_core, collect(counts))
    ncore = length(ix_core) * length(iy_core)
    occ = ncore == 0 ? 0.0 : incore / ncore
    return (maxper = maxper, occupied = occ, outside = outside)
end

#----- local frame to the DEM's coordinate system -----
function _data_origin_latlon(d::Data)
    ok = isfinite.(d.lat) .& isfinite.(d.lon)
    any(ok) || error("the stations have no latitude/longitude (frame coordinates): give the ground as `surface` instead of a topography file")
    lat0 = (length(d.origin) >= 1 && isfinite(d.origin[1])) ? Float64(d.origin[1]) : mean(d.lat[ok])
    lon0 = (length(d.origin) >= 2 && isfinite(d.origin[2])) ? Float64(d.origin[2]) : mean(d.lon[ok])
    return lat0, lon0
end

function _local_tm_to_target_context(d::Data, target_crs::AbstractString)
    _data_origin_latlon(d)                                       # stations must carry latitudes/longitudes
    # a TM centred on the stations, fed frame offsets from their centroid: whatever the frame (a local TM, a
    # national grid), it then differs from this TM by a near-rigid rotation, scale and shift, taken up below
    ok = isfinite.(d.lat) .& isfinite.(d.lon)
    lat0, lon0 = mean(d.lat[ok]), mean(d.lon[ok])
    fx, fy = mean(d.x[ok]), mean(d.y[ok])
    src = "+proj=tmerc +lat_0=$(float(lat0)) +lon_0=$(float(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs"
    local_tm_to_target = Proj.Transformation(src, strip(target_crs); always_xy = true)
    wgs84_to_target = Proj.Transformation("EPSG:4326", strip(target_crs); always_xy = true)

    pred_x = Float64[]
    pred_y = Float64[]
    actual_x = Float64[]
    actual_y = Float64[]
    for i in eachindex(d.x)
        isfinite(d.lat[i]) && isfinite(d.lon[i]) || continue
        p = local_tm_to_target((Float64(d.y[i]) - fy + 500000.0, Float64(d.x[i]) - fx))
        a = wgs84_to_target((Float64(d.lon[i]), Float64(d.lat[i])))
        push!(pred_x, Float64(p[1]))
        push!(pred_y, Float64(p[2]))
        push!(actual_x, Float64(a[1]))
        push!(actual_y, Float64(a[2]))
    end

    # an affine map, least squares over the stations, takes the TM prediction to the stations' true target
    # coordinates (rotation, scale and shift). Centred on the stations for conditioning
    cx, cy = mean(pred_x), mean(pred_y)
    P = hcat(pred_x .- cx, pred_y .- cy, ones(length(pred_x)))
    M = length(pred_x) >= 3 && rank(P) == 3 ? permutedims(P \ hcat(actual_x, actual_y)) :
        [1.0 0.0 median(actual_x .- pred_x); 0.0 1.0 median(actual_y .- pred_y)]
    return (
        trans = local_tm_to_target,
        affine = M, centre = (cx, cy), frame_centre = (fx, fy),
        lat0 = lat0,
        lon0 = lon0,
    )
end

function _local_xy_to_target_xy(ctx, x_local::Real, y_local::Real)
    p = ctx.trans((Float64(y_local) - ctx.frame_centre[2] + 500000.0, Float64(x_local) - ctx.frame_centre[1]))
    t = ctx.affine * [p[1] - ctx.centre[1], p[2] - ctx.centre[2], 1.0]
    return t[1], t[2]
end

# the inverse: target CRS coordinates to frame (x north, y east)
function _target_xy_to_local_xy(ctx, back, tx::Real, ty::Real)
    q = ctx.affine[:, 1:2] \ ([tx, ty] .- ctx.affine[:, 3])
    p = back((q[1] + ctx.centre[1], q[2] + ctx.centre[2]))
    return Float64(p[2]) + ctx.frame_centre[1], Float64(p[1]) - 500000.0 + ctx.frame_centre[2]
end

function _load_topography_geotiff(path::AbstractString)
    ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        # ArchGDAL returns (x, y); sampling below indexes zgrid[y, x].
        zgrid = permutedims(Float64.(ArchGDAL.read(band)))
        nodata = ArchGDAL.getnodatavalue(band)
        gt = ArchGDAL.getgeotransform(dataset)
        x0, dx, rx, y0, ry, dy = Float64.(gt)

        (abs(rx) > 1e-9 || abs(ry) > 1e-9) && error("Rotated GeoTIFF geotransforms are not supported by this mesh builder.")

        ny, nx = size(zgrid)
        xs = x0 .+ ((0:(nx - 1)) .+ 0.5) .* dx
        ys = y0 .+ ((0:(ny - 1)) .+ 0.5) .* dy

        if xs[end] < xs[1]
            xs = reverse(xs)
            zgrid = reverse(zgrid; dims = 2)
        end
        if ys[end] < ys[1]
            ys = reverse(ys)
            zgrid = reverse(zgrid; dims = 1)
        end

        if nodata !== nothing
            nd = Float64(nodata)
            zgrid[zgrid .== nd] .= NaN
        end
        # Some DEM products use extreme finite sentinels instead of explicit NaN.
        zgrid[abs.(zgrid) .> 1e30] .= NaN
        zgrid[.!isfinite.(zgrid)] .= NaN

        valid_rows = vec(any(isfinite.(zgrid), dims = 2))
        valid_cols = vec(any(isfinite.(zgrid), dims = 1))
        any(valid_rows) || error("Topography GeoTIFF contains no finite elevation samples.")
        any(valid_cols) || error("Topography GeoTIFF contains no finite elevation samples.")
        valid_bbox = (
            minimum(xs[valid_cols]), maximum(xs[valid_cols]),
            minimum(ys[valid_rows]), maximum(ys[valid_rows]),
        )

        return (
            x = collect(xs),
            y = collect(ys),
            z = zgrid,
            bbox = (minimum(xs), maximum(xs), minimum(ys), maximum(ys)),
            valid_bbox = valid_bbox,
        )
    end
end

function _ensure_topography(path::AbstractString)
    haskey(_TOPO_CACHE, path) && return _TOPO_CACHE[path]
    isfile(path) || error("Topography file not found: $path")
    return _TOPO_CACHE[path] = _load_topography_geotiff(path)
end

# inside the DEM's valid samples (the sampler below clamps outside points to the edge, for padding columns only)
_in_dem(topo, x, y) = (b = topo.valid_bbox; b[1] <= x <= b[2] && b[3] <= y <= b[4])

function _sample_topography(topo, xq::Real, yq::Real)
    xs = topo.x
    ys = topo.y
    z = topo.z
    nx = length(xs)
    ny = length(ys)

    # Padding cells can extend beyond the DEM. Clamp them to the raster edge so
    # the topography over the padding decays to the nearest available elevation.
    xlo, xhi, ylo, yhi = topo.valid_bbox
    xqc = clamp(Float64(xq), xlo, xhi)
    yqc = clamp(Float64(yq), ylo, yhi)

    ix = xqc == xs[end] ? nx - 1 : clamp(searchsortedlast(xs, xqc), 1, nx - 1)
    iy = yqc == ys[end] ? ny - 1 : clamp(searchsortedlast(ys, yqc), 1, ny - 1)

    x1, x2 = xs[ix], xs[ix + 1]
    y1, y2 = ys[iy], ys[iy + 1]
    tx = x2 == x1 ? 0.0 : (xqc - x1) / (x2 - x1)
    ty = y2 == y1 ? 0.0 : (yqc - y1) / (y2 - y1)

    z11 = z[iy, ix]
    z12 = z[iy, ix + 1]
    z21 = z[iy + 1, ix]
    z22 = z[iy + 1, ix + 1]
    if all(isfinite, (z11, z12, z21, z22))
        return (1 - tx) * (1 - ty) * z11 + tx * (1 - ty) * z12 + (1 - tx) * ty * z21 + tx * ty * z22
    end

    candidates = [(z11, hypot(xqc - x1, yqc - y1)), (z12, hypot(xqc - x2, yqc - y1)),
                  (z21, hypot(xqc - x1, yqc - y2)), (z22, hypot(xqc - x2, yqc - y2))]
    finite_candidates = filter(c -> isfinite(c[1]), candidates)
    if !isempty(finite_candidates)
        return first(finite_candidates[argmin(last.(finite_candidates))])
    end

    # If the immediate bilinear stencil falls in a no-data moat, expand locally
    # until the nearest finite DEM sample is found.
    for radius in 1:25
        ix1 = max(1, ix - radius)
        ix2 = min(nx, ix + radius + 1)
        iy1 = max(1, iy - radius)
        iy2 = min(ny, iy + radius + 1)
        best_val = NaN
        best_dist = Inf
        for jy in iy1:iy2, jx in ix1:ix2
            val = z[jy, jx]
            isfinite(val) || continue
            dist = hypot(xqc - xs[jx], yqc - ys[jy])
            if dist < best_dist
                best_dist = dist
                best_val = val
            end
        end
        isfinite(best_val) && return best_val
    end
    return NaN
end

# ground depth (z down, sea level 0, the frame of the station heights) at every column centre, from the
# DEM (elevations in m above the datum of the station heights) or from a function of (x, y)
function _sample_surface_local(m, d, ctx, topo_path::AbstractString)
    topo = _ensure_topography(topo_path)
    station_elev = Float64[]
    outside = count(i -> !_in_dem(topo, _local_xy_to_target_xy(ctx, d.x[i], d.y[i])...), eachindex(d.x))
    outside == 0 || error("$outside of $(length(d.x)) stations lie outside the topography file $(basename(topo_path))")
    for i in eachindex(d.x)
        tx, ty = _local_xy_to_target_xy(ctx, d.x[i], d.y[i])
        e = _sample_topography(topo, tx, ty)
        isfinite(e) && push!(station_elev, e)
    end
    isempty(station_elev) && error("The topography file does not cover the stations.")
    fallback_elev = median(station_elev)

    surface_local = Matrix{Float64}(undef, m.nx, m.ny)
    fallback_count = 0
    for ix in 1:m.nx, iy in 1:m.ny
        tx, ty = _local_xy_to_target_xy(ctx, m.x_centers[ix], m.y_centers[iy])
        elev_abs = _sample_topography(topo, tx, ty)
        if !isfinite(elev_abs)
            elev_abs = fallback_elev
            fallback_count += 1
        end
        surface_local[ix, iy] = -elev_abs
    end

    fallback_count == 0 || @warn("Used median DEM elevation fallback for padded mesh columns that still sampled nodata.", fallback_count = fallback_count)

    return surface_local
end

_surface_from_function(m, surface) = [-Float64(surface(x, y)) for x in m.x_centers, y in m.y_centers]

# padding widths growing by `factor` from `h` until `extent` is covered, or `n` cells
function _padding_widths(h, factor, extent, n)
    n !== nothing && return [h * factor^i for i in 1:Int(n)]
    p = Float64[]
    while sum(p; init = 0.0) < extent
        push!(p, h * factor^(length(p) + 1))
    end
    p
end

#----- build the mesh from stations, padding and depth settings -----
function design_mesh(sx, sy; ρ_bg = 0.0, dx_core, dy_core, nx_pad = nothing, ny_pad = nothing,
                     pad_factor, z_first, z_factor, depth_mult, uniform_layers = nothing,
                     nx_core = nothing, ny_core = nothing,
                     surface_core_range = nothing, n_air = 0,
                     air_factor = 1.35, air_first = nothing,
                     air_first_div = 2.0)
    dxc = max(1e-3, dx_core)
    dyc = max(1e-3, dy_core)

    # one cell beyond the stations; an odd count puts the survey centre (and a regular grid of
    # stations at twice the cell size) on cell centres, off the mesh planes
    odd(n) = n + iseven(n)
    cx_mid = (minimum(sx) + maximum(sx)) / 2
    cy_mid = (minimum(sy) + maximum(sy)) / 2
    nx_core = nx_core === nothing ? odd(ceil(Int, (maximum(sx) - minimum(sx)) / dxc) + 2) : max(1, Int(nx_core))
    ny_core = ny_core === nothing ? odd(ceil(Int, (maximum(sy) - minimum(sy)) / dyc) + 2) : max(1, Int(ny_core))
    x0 = cx_mid - nx_core * dxc / 2
    y0 = cy_mid - ny_core * dyc / 2

    padx = _padding_widths(dxc, pad_factor, nx_core * dxc, nx_pad === nothing ? nothing : max(0, nx_pad))
    pady = _padding_widths(dyc, pad_factor, ny_core * dyc, ny_pad === nothing ? nothing : max(0, ny_pad))
    nx_pad, ny_pad = length(padx), length(pady)
    dx = vcat(reverse(padx), fill(dxc, nx_core), padx)
    dy = vcat(reverse(pady), fill(dyc, ny_core), pady)

    origin_x = x0 - sum(padx; init = 0.0)
    origin_y = y0 - sum(pady; init = 0.0)

    # regular layers from the highest ground of the core to the depth of investigation below its lowest
    # ground, then layers growing by z_factor until they are as thick as the regular block. With
    # `uniform_layers` = N the regular block ends N layers below the lowest ground (lake beds included), and
    # the growing layers carry on to the same total depth, so thin layers can follow the relief and the lakes
    doi = depth_mult * min(nx_core * dxc, ny_core * dyc)
    top, bottom = surface_core_range === nothing ? (0.0, 0.0) : Float64.(surface_core_range)
    if uniform_layers === nothing
        block = fill(Float64(z_first), max(1, ceil(Int, (bottom - top + doi) / z_first - 1e-9)))
        ground_dz = vcat(block, _padding_widths(Float64(z_first), z_factor, sum(block), nothing))
    else
        Int(uniform_layers) >= 0 || throw(ArgumentError("uniform_layers must be nonnegative"))
        block = fill(Float64(z_first), max(1, ceil(Int, (bottom - top) / z_first - 1e-9) + Int(uniform_layers)))
        ground_dz = vcat(block, _padding_widths(Float64(z_first), z_factor, 2 * (bottom - top + doi) - sum(block), nothing))
    end

    air_dz = Float64[]
    if Int(n_air) > 0
        air_seed = air_first === nothing ? max(1.0, z_first / air_first_div) : max(1.0, Float64(air_first))
        air_dz = reverse([air_seed * air_factor^(i - 1) for i in 1:Int(n_air)])
    end
    origin_z = top - sum(air_dz; init = 0.0)

    dz = vcat(air_dz, ground_dz)
    x_edges = origin_x .+ vcat(0.0, cumsum(dx))
    y_edges = origin_y .+ vcat(0.0, cumsum(dy))
    z_edges = origin_z .+ vcat(0.0, cumsum(dz))
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2
    z_centers = (z_edges[1:end-1] .+ z_edges[2:end]) ./ 2

    ix_core = (nx_pad + 1):(nx_pad + nx_core)
    iy_core = (ny_pad + 1):(ny_pad + ny_core)
    diag = site_occupancy(sx, sy, x_edges, y_edges, ix_core, iy_core)

    core_x0 = x_edges[nx_pad + 1]
    core_x1 = x_edges[nx_pad + nx_core + 1]
    core_y0 = y_edges[ny_pad + 1]
    core_y1 = y_edges[ny_pad + ny_core + 1]
    pad_extent = min(sum(padx; init = 0.0), sum(pady; init = 0.0))
    core_width = min(core_x1 - core_x0, core_y1 - core_y0)
    nz_air = length(air_dz)
    nz_earth = length(ground_dz)

    return (
        dx = dx, dy = dy, dz = dz, origin = [origin_x, origin_y, origin_z], ρ_bg = Float64(ρ_bg),
        nx = length(dx), ny = length(dy), nz = length(dz), nz_air = nz_air, nz_earth = nz_earth, nz_regular = length(block),
        nx_core = nx_core, ny_core = ny_core, nx_pad = nx_pad, ny_pad = ny_pad, dx_core = dxc, dy_core = dyc,
        x_edges = x_edges, y_edges = y_edges, z_edges = z_edges,
        x_centers = x_centers, y_centers = y_centers, z_centers = z_centers,
        x_edges_km = x_edges ./ 1000, y_edges_km = y_edges ./ 1000, z_edges_km = z_edges ./ 1000,
        core_x0_km = core_x0 / 1000, core_x1_km = core_x1 / 1000,
        core_y0_km = core_y0 / 1000, core_y1_km = core_y1 / 1000,
        full_x0_km = x_edges[1] / 1000, full_x1_km = x_edges[end] / 1000,
        full_y0_km = y_edges[1] / 1000, full_y1_km = y_edges[end] / 1000,
        doi_km = doi / 1000, core_width_km = core_width / 1000,
        depth_km = sum(ground_dz) / 1000, pad_extent_km = pad_extent / 1000,
        surface_min_km = top / 1000, surface_max_km = bottom / 1000,
        pad_ok = pad_extent >= core_width,
        maxper = diag.maxper, occupied = diag.occupied, outside = diag.outside,
    )
end

function _build_density_model(m, surface_local)
    A = fill(Float64(m.ρ_bg), m.nx, m.ny, m.nz)
    @inbounds for ix in 1:m.nx, iy in 1:m.ny, iz in 1:m.nz
        ground = surface_local === nothing ? m.origin[3] + sum(m.dz[1:m.nz_air]) : surface_local[ix, iy]
        if m.z_centers[iz] < ground
            A[ix, iy, iz] = NaN
        end
    end
    return A
end

# station column: the deepest station sets how far down the ground must start
function _site_column_constraints(m, d)
    cols = Dict{Tuple{Int, Int}, Tuple{String, Float64}}()
    for isite in eachindex(d.x)
        ix = clamp(searchsortedlast(m.x_edges, d.x[isite]), 1, m.nx)
        iy = clamp(searchsortedlast(m.y_edges, d.y[isite]), 1, m.ny)
        key = (ix, iy)
        site_z = Float64(d.z[isite])
        current = get(cols, key, ("", -Inf))
        if site_z > current[2]
            cols[key] = (String(d.code[isite]), site_z)
        end
    end
    return cols
end

# lower the ground in station columns until the top of the first earth cell lies at or below each station
function _apply_site_surface_constraints(m, d, surface_local)
    adjusted = copy(surface_local)
    adjustments = NamedTuple[]
    impossible = NamedTuple[]

    for ((ix, iy), (site_name, site_z)) in _site_column_constraints(m, d)
        k = findfirst(>=(site_z), m.z_edges[1:end-1])        # first cell whose top lies at or below the station
        if k === nothing
            push!(impossible, (site = site_name, ix = ix, iy = iy, site_z = site_z,
                               bottom = Float64(m.z_edges[end]), clearance = Float64(m.z_edges[end]) - site_z))
            continue
        end
        old_surface = Float64(adjusted[ix, iy])
        first_earth = something(findfirst(>=(old_surface), m.z_centers), m.nz + 1)
        first_earth >= k && continue                          # the station already lies above the ground
        new_surface = m.z_centers[k]
        adjusted[ix, iy] = new_surface
        if new_surface > old_surface + 1e-9
            push!(adjustments, (site = site_name, ix = ix, iy = iy, site_z = site_z,
                                surface_from = old_surface, surface_to = new_surface, drop = new_surface - old_surface))
        end
    end

    sort!(adjustments; by = row -> -row.drop)
    sort!(impossible; by = row -> row.clearance)
    text = isempty(adjustments) ?
        "The ground follows the topography only." :
        @sprintf("Ground lowered in %d station column(s) so every station lies on or above it.", length(adjustments))

    return adjusted, (nadjusted = length(adjustments), adjustments = adjustments, impossible = impossible, text = text)
end

function _site_ground_report(m, d)
    bad = NamedTuple[]
    for isite in eachindex(d.x)
        ix = clamp(searchsortedlast(m.x_edges, d.x[isite]), 1, m.nx)
        iy = clamp(searchsortedlast(m.y_edges, d.y[isite]), 1, m.ny)
        first_earth = findfirst(!isnan, vec(m.A[ix, iy, :]))
        first_earth === nothing && error("No earth cells remain in the column of station $(d.code[isite]) at mesh column ($ix, $iy).")
        top_edge = m.z_edges[first_earth]
        clearance = Float64(top_edge) - Float64(d.z[isite])
        if clearance < 0.0
            push!(bad, (site = String(d.code[isite]), ix = ix, iy = iy, site_z = Float64(d.z[isite]),
                        top_edge = Float64(top_edge), clearance = clearance))
        end
    end

    sort!(bad; by = row -> row.clearance)
    nbad = length(bad)
    ok = nbad == 0
    sample_rows = bad[1:min(end, 5)]
    sample = isempty(sample_rows) ? "" : join([@sprintf("%s needs %.1f m", row.site, -row.clearance) for row in sample_rows], ", ")
    max_drop = isempty(bad) ? 0.0 : -bad[1].clearance
    summary = ok ?
        @sprintf("✓ All %d stations lie on or above the first earth cell.", length(d.x)) :
        @sprintf("• %d station(s) lie below the top of the first earth cell. Example: %s.", nbad, sample)
    fix = ok ? "" :
        "Fix the mesh, not the data: lower the ground (or thin the first layer) in the flagged columns until every station is on or above the top of its first earth cell."

    return (ok = ok, nbad = nbad, rows = bad, sample = sample, max_drop = max_drop, text = summary, fix = fix)
end

function _print_site_surface_report(m)
    if m.site_surface.nadjusted > 0
        println(m.site_surface.text)
        for row in m.site_surface.adjustments[1:min(end, 10)]
            println(@sprintf("  %s: xcell=%d ycell=%d ground %.3f -> %.3f (lowered %.3f m)",
                             row.site, row.ix, row.iy, row.surface_from, row.surface_to, row.drop))
        end
        if m.site_surface.nadjusted > 10
            println(@sprintf("  ... %d more adjusted column(s) omitted.", m.site_surface.nadjusted - 10))
        end
    end

    if !isempty(m.site_surface.impossible)
        @warn(@sprintf("%d station(s) lie below the bottom of the mesh; revise the vertical mesh or the station heights.", length(m.site_surface.impossible)))
        for row in m.site_surface.impossible[1:min(end, 10)]
            println(@sprintf("  %s: xcell=%d ycell=%d site_z=%.3f bottom=%.3f", row.site, row.ix, row.iy, row.site_z, row.bottom))
        end
    end
    return nothing
end

function _print_site_ground_report(m)
    m.site_ground.ok && return nothing
    @warn(@sprintf("%d station(s) lie at or below the ground of the written model.", m.site_ground.nbad))
    println(m.site_ground.fix)
    for row in m.site_ground.rows[1:min(end, 10)]
        println(@sprintf("  %s: xcell=%d ycell=%d site_z=%.3f top_of_first_earth=%.3f clearance=%.3f",
                         row.site, row.ix, row.iy, row.site_z, row.top_edge, row.clearance))
    end
    m.site_ground.nbad > 10 && println(@sprintf("  ... %d more station(s) omitted.", m.site_ground.nbad - 10))
    return nothing
end

function build_mesh_bundle(d, sx, sy;
                           ρ_bg = 0.0, ρ_ref = 2670.0, dx_core, dy_core, nx_core, ny_core,
                           nx_pad = nothing, ny_pad = nothing, pad_factor, z_first, z_factor, depth_mult,
                           uniform_layers = nothing, topo_ctx = nothing, topo_path::AbstractString = "", surface = nothing,
                           lakes = NamedTuple[], n_air::Int = 0, air_factor::Real = 1.35, air_first_div::Real = 2.0)
    use_topo = !isempty(strip(topo_path))
    use_topo && surface !== nothing && error("give the ground as a topography file or as `surface`, not both")
    kw = (ρ_bg = Float64(ρ_bg), dx_core = Float64(dx_core), dy_core = Float64(dy_core), nx_core = Int(nx_core),
          ny_core = Int(ny_core), nx_pad = nx_pad, ny_pad = ny_pad, pad_factor = Float64(pad_factor),
          z_first = Float64(z_first), z_factor = Float64(z_factor), depth_mult = Float64(depth_mult),
          uniform_layers = uniform_layers)
    base = design_mesh(sx, sy; kw...)

    # lakes (in the frame, see Lakes3D.jl) lay their beds into the ground before the layers are fitted to it
    surface_local = nothing
    water_level = fill(NaN, base.nx, base.ny)
    if use_topo || surface !== nothing || !isempty(lakes)
        use_topo && topo_ctx === nothing && error("Topography mode was requested without a valid projection context.")
        flat(b) = fill(0.0, b.nx, b.ny)
        sample(b) = use_topo ? _sample_surface_local(b, d, topo_ctx, topo_path) :
                    surface !== nothing ? _surface_from_function(b, surface) : flat(b)
        water_level, s0 = _apply_lakes(base.x_centers, base.y_centers, sample(base), lakes)
        core = s0[base.nx_pad+1:base.nx_pad+base.nx_core, base.ny_pad+1:base.ny_pad+base.ny_core]
        # the lateral mesh does not depend on the layers, so the ground sampled on it holds for the final mesh
        base = design_mesh(sx, sy; kw..., surface_core_range = extrema(core), n_air = n_air,
                           air_factor = Float64(air_factor), air_first_div = Float64(air_first_div))
        surface_local = s0
    elseif n_air > 0
        base = design_mesh(sx, sy; kw..., n_air = n_air, air_factor = Float64(air_factor), air_first_div = Float64(air_first_div))
    end

    site_surface = (nadjusted = 0, adjustments = NamedTuple[], impossible = NamedTuple[],
                    text = "The ground follows the topography only.")
    if surface_local === nothing
        surface_local = fill(base.origin[3] + sum(base.dz[1:base.nz_air]; init = 0.0), base.nx, base.ny)
    end
    _check_stations_above_water(base, d, water_level)
    surface_local, site_surface = _apply_site_surface_constraints(base, d, surface_local)

    A = _build_density_model(base, surface_local)
    # water: the air cells of a lake column between its level and its bed, at the water contrast 1000 − ρ_ref, fixed
    ρ_water = 1000.0 - Float64(ρ_ref)
    water = falses(size(A))
    for ix in 1:base.nx, iy in 1:base.ny
        isnan(water_level[ix, iy]) && continue
        for iz in 1:base.nz
            isnan(A[ix, iy, iz]) && base.z_centers[iz] > water_level[ix, iy] && (A[ix, iy, iz] = ρ_water; water[ix, iy, iz] = true)
        end
    end
    site_ground = _site_ground_report((; base..., A = A), d)
    # the mesh is laid out in the frame of the data, so it carries the data's azimuth as its rotation
    return (; base..., A = A, mask = .!isnan.(A) .& .!water, water = water, nlakes = length(lakes), rotation = d.azimuth,
              topo_active = use_topo || surface !== nothing,
              topo_name = use_topo ? basename(topo_path) : surface !== nothing ? "surface function" : "flat",
              surface_local = surface_local,
              site_surface = site_surface,
              site_ground = site_ground)
end

"""
    mesh_model(m) -> Model

The start model of a mesh bundle (`MakeMesh3D`, `build_mesh_bundle`) as a `Model`, with z = 0 at the highest
ground: the top of the highest earth cell, whose elevation becomes the model's `datum`, and the rotation of the
data frame (the data's azimuth). The bundle itself is in the frame of the data (z = −elevation).
"""
function mesh_model(m)
    k = something(findfirst(k -> any(!isnan, view(m.A, :, :, k)), 1:m.nz), 1)
    ztop = m.origin[3] + sum(m.dz[1:k-1]; init = 0.0)                 # top of the highest earth cell, z = −elevation
    Model(m.dx, m.dy, m.dz, m.A; origin = [m.origin[1], m.origin[2], m.origin[3] - ztop], rotation = m.rotation, datum = -ztop)
end

# no station may lie below the level of the lake under it (a station over water flies above it)
function _check_stations_above_water(m, d, level)
    for n in eachindex(d.x)
        ix = clamp(searchsortedlast(m.x_edges, d.x[n]), 1, m.nx)
        iy = clamp(searchsortedlast(m.y_edges, d.y[n]), 1, m.ny)
        z = level[ix, iy]
        isnan(z) || d.z[n] <= z ||
            error(@sprintf("station %s lies %.1f m below the level of the lake under it", d.code[n], d.z[n] - z))
    end
end

function _save_mesh_outputs(m, out_model::AbstractString, out_mask::AbstractString)
    model = mesh_model(m)
    write_model(out_model, model; comment = "start model from MakeMesh3D")
    write_mask(out_mask, model, m.mask; water = m.water)
    _print_site_surface_report(m)
    _print_site_ground_report(m)
    return nothing
end

"""
    MakeMesh3D(data_file; out_model, out_mask, topo_file="", topo_crs="EPSG:4326", surface=nothing,
               lakes=[], lakes_crs="EPSG:4326", mode=:gui, cell_width_frac=0.5, n_pad=nothing, pad_factor=1.4,
               first_layer_div=2.0, vertical_factor=1.4, depth_mult=0.5, uniform_layers=nothing, ρ_bg=0.0,
               ρ_ref=2670.0, air_layers=0, ...)

Design a mesh for a gravity data file (native format, stations in the mesh frame) and write the start
model (`mesh_start_model.den`, density contrast `ρ_bg`, air NaN) and the inversion mask
(`mesh_start_model.mask`) next to the data.

- core cells `cell_width_frac` × the median nearest-station spacing, an odd count one cell beyond the stations
- `n_pad` padding cells per side growing by `pad_factor` (nothing: until the padding is as wide as the core)
- regular layers of spacing / `first_layer_div` from the highest ground of the core to the depth of
  investigation `depth_mult` × core width below its lowest ground, then layers growing by `vertical_factor`
  until they are as thick as the regular block; with `uniform_layers` = N the regular layers end N layers below
  the lowest ground (lake beds included) and grow from there to the same total depth (thin layers for the relief
  and the lakes at a fraction of the cells)
- the ground from `topo_file` (GeoTIFF, elevations in `topo_crs`, needs station latitudes/longitudes) or
  `surface(x, y)` (elevation, m; frame coordinates); flat at the top of the earth otherwise. Columns are
  lowered where a station would lie below the ground
- `lakes`: `[(polygon, level[, depth, taper]), …]` in `lakes_crs` ("frame" for mesh coordinates; see `Lakes3D.jl`,
  `read_lakes`): cells between a
  lake's level and its bed are water, the contrast 1000 − `ρ_ref` (kg/m³; `ρ_ref` the reduction density of the
  Bouguer anomaly) and mask 9 (fixed, outside the smoothing); the ground is the bed
  unless the lake gives a `depth`. Stations may lie over a lake, not below its level
- `air_layers` air cells above the highest ground, thinnest at the ground

`mode = :gui` opens the designer window (GLMakie; Save writes the files) and returns `(fig, mesh)`;
`interactive = false` only builds it. `:nogui` writes at once and returns the mesh bundle (`mesh_model`
turns it into a `Model`).
"""
function MakeMesh3D(data_file::AbstractString;
    out_model::AbstractString = joinpath(dirname(abspath(data_file)), "mesh_start_model.den"),
    out_mask::AbstractString = joinpath(dirname(out_model), "mesh_start_model.mask"),
    topo_file::AbstractString = "",
    topo_crs::AbstractString = "EPSG:4326",
    surface = nothing,
    lakes = NamedTuple[],
    lakes_crs::AbstractString = "EPSG:4326",
    ρ_bg::Real = 0.0,
    ρ_ref::Real = 2670.0,
    mode::Symbol = :gui,
    cell_width_frac::Real = 0.5,
    n_pad::Union{Nothing, Int} = nothing,
    pad_factor::Real = 1.4,
    first_layer_div::Real = 2.0,
    vertical_factor::Real = 1.4,
    depth_mult::Real = 0.5,
    uniform_layers::Union{Nothing, Int} = nothing,
    dx_core = nothing,
    dy_core = nothing,
    z_first = nothing,
    air_layers::Int = 0,
    air_factor::Real = 1.35,
    air_first_div::Real = 2.0,
    colormap = COLORMAP,
    density_range = (-300.0, 300.0),
    site_color = :black,
    site_size_full = 6,
    site_size_core = 5,
    grid_color = (:grey, 0.7),
    grid_linewidth = 1.0,
    fig_size = (1500, 1000),
    interactive::Bool = true)

    mode in (:gui, :nogui) || throw(ArgumentError("mode must be :gui or :nogui, got $mode"))
    isempty(data_file) && error("data_file is required")
    isfile(data_file) || error("Data file not found: $data_file")

    use_topo = !isempty(strip(topo_file))
    use_topo && (isfile(topo_file) || error("Topography file not found: $topo_file"))

    d = load_data_grav(data_file)
    topo_ctx = use_topo ? _local_tm_to_target_context(d, topo_crs) : nothing
    lakes_frame = _lakes_to_frame(d, lakes, lakes_crs)
    sx = collect(Float64.(d.x))
    sy = collect(Float64.(d.y))

    spacing = nearest_neighbour_spacing(sx, sy)
    nice(v) = round(v; sigdigits = 3)
    dx_core0 = Float64(something(dx_core, nice(spacing.median * cell_width_frac)))
    dy_core0 = Float64(something(dy_core, dx_core0))
    odd(n) = n + iseven(n)
    nx_core0 = odd(ceil(Int, (maximum(sx) - minimum(sx)) / dx_core0) + 2)
    ny_core0 = odd(ceil(Int, (maximum(sy) - minimum(sy)) / dy_core0) + 2)
    z_first0 = Float64(something(z_first, nice(spacing.median / first_layer_div)))

    @info @sprintf("Loaded %d stations (%s); median spacing %.0f m; core cell %.0f m, layer %.0f m; topography %s; %d lake(s)",
                   length(sx), _COMPONENT_NAME[d.component], spacing.median, dx_core0, z_first0,
                   use_topo ? basename(topo_file) : surface !== nothing ? "surface function" : "off", length(lakes_frame))

    if mode === :nogui
        m = build_mesh_bundle(d, sx, sy;
            ρ_bg = Float64(ρ_bg), ρ_ref = Float64(ρ_ref), dx_core = dx_core0, dy_core = dy_core0,
            nx_core = nx_core0, ny_core = ny_core0, nx_pad = n_pad, ny_pad = n_pad,
            pad_factor = pad_factor, z_first = z_first0, z_factor = vertical_factor,
            depth_mult = depth_mult, uniform_layers = uniform_layers, topo_ctx = topo_ctx, topo_path = topo_file,
            surface = surface, lakes = lakes_frame, n_air = air_layers, air_factor = air_factor, air_first_div = air_first_div)
        @printf("grid %d×%d×%d (%d cells); core %d×%d @ %.0f×%.0f m; air %d; %d regular layers; water %d cells; max stations/cell %d; occupied %.0f%%; pad %.2f km/side; depth %.2f km (investigation %.2f km)\n",
            m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core,
            m.nz_air, m.nz_regular, count(m.water), m.maxper, 100 * m.occupied, m.pad_extent_km, m.depth_km, m.doi_km)
        foreach(s -> println(s), suggestions(m).lines)
        _save_mesh_outputs(m, out_model, out_mask)
        @printf("wrote %s\n", out_model)
        @printf("wrote %s\n", out_mask)
        return m
    end

    isdefined(GravityGeophysics, :_make_mesh3D_gui) ||
        error("GLMakie is not available; use mode = :nogui on headless systems")

    return _make_mesh3D_gui(;
        data_file = data_file,
        out_model = out_model,
        out_mask = out_mask,
        topo_file = topo_file,
        surface = surface,
        lakes = lakes_frame,
        uniform_layers = uniform_layers,
        n_pad = n_pad,
        pad_factor = pad_factor,
        vertical_factor = vertical_factor,
        depth_mult = depth_mult,
        air_layers = air_layers,
        air_factor = air_factor,
        air_first_div = air_first_div,
        colormap = colormap,
        density_range = density_range,
        site_color = site_color,
        site_size_full = site_size_full,
        site_size_core = site_size_core,
        grid_color = grid_color,
        grid_linewidth = grid_linewidth,
        fig_size = fig_size,
        interactive = interactive,
        d = d, sx = sx, sy = sy,
        spacing = spacing, ρ_bg0 = Float64(ρ_bg), ρ_ref = Float64(ρ_ref),
        dx_core0 = dx_core0, dy_core0 = dy_core0,
        nx_core0 = nx_core0, ny_core0 = ny_core0, z_first0 = z_first0,
        topo_ctx = topo_ctx)
end

#----- check the mesh and advise what to fix -----
function suggestions(m)
    msgs = String[]
    m.maxper > 1 && push!(msgs,
        "• More than one station falls in some cells — decrease 'Core cell N–S (m)' / 'Core cell E–W (m)'.")
    m.outside > 0 && push!(msgs,
        "• $(m.outside) station(s) lie outside the core — increase 'Core cells N–S' / 'Core cells E–W'.")
    m.pad_ok || push!(msgs,
        "• Lateral padding is narrower than the core — increase 'Pad cells N–S' / 'Pad cells E–W' or 'Pad factor'.")
    !isempty(m.site_surface.impossible) && push!(msgs,
        "• Some stations lie below the bottom of the mesh — check the station heights and the topography datum.")
    m.site_ground.ok || push!(msgs, m.site_ground.text * " " * m.site_ground.fix)
    ok = isempty(msgs)
    text = ok ? "✓ Mesh looks well-sized for gravity inversion." : join(msgs, "    ")
    return (text = text, ok = ok, lines = msgs)
end
