# 3D mesh design for gravity inversion
# Author: @pankajkmishra
# Tensor mesh from the survey: uniform core under the data, geometric lateral padding, layers growing
# with depth, optional topography (cells above the ground are air, covariance mask 0) and the covariance
# Depth is measured down from the highest ground, so the data carry an elevation rather than a negative z,
# and the topography is a callable or a flat elevation, which keeps a GeoTIFF reader out of the dependencies

using Printf
using Statistics

_snap(v, r) = collect(r)[argmin(abs.(collect(r) .- v))]

#----- find the typical data spacing -----

"""
    data_spacing(sx, sy) -> (median, mean)

Typical station spacing: the median nearest-neighbour distance and `sqrt(area / n)`.
"""
function data_spacing(sx, sy)
    n = length(sx)
    n < 2 && return (median = 100.0, mean = 100.0)
    area = max(maximum(sx) - minimum(sx), 1.0) * max(maximum(sy) - minimum(sy), 1.0)
    # nearest neighbour on a coarse bucket grid, so dense surveys stay cheap
    h       = sqrt(area / n)
    buckets = Dict{Tuple{Int, Int}, Vector{Int}}()
    for i in 1:n
        push!(get!(buckets, (floor(Int, sx[i] / h), floor(Int, sy[i] / h)), Int[]), i)
    end
    dmin = fill(Inf, n)
    for i in 1:n
        bx, by = floor(Int, sx[i] / h), floor(Int, sy[i] / h)
        for ox in -1:1, oy in -1:1, j in get(buckets, (bx + ox, by + oy), Int[])
            j == i || (dmin[i] = min(dmin[i], hypot(sx[i] - sx[j], sy[i] - sy[j])))
        end
    end
    d = filter(isfinite, dmin)
    (median = isempty(d) ? h : median(d), mean = h)
end

#----- build the mesh from data, padding and depth settings -----

"""
    design_mesh(sx, sy; rho_bg, dx_core, dy_core, nx_pad, ny_pad, pad_factor, z_first, z_factor,
                depth_mult, nx_core=nothing, ny_core=nothing, relief=0.0)

Tensor mesh around the data at north `sx` and east `sy`: a uniform core centred on the data, `n_pad` padding cells growing by `pad_factor` on each side, and layers from the highest ground down, `z_first` thick through the relief (`relief`, the drop from the highest to the lowest ground), then growing by `z_factor` to `depth_mult` times the narrower core width below the lowest ground.
"""
function design_mesh(sx, sy; rho_bg, dx_core, dy_core, nx_pad, ny_pad, pad_factor, z_first, z_factor,
                     depth_mult, nx_core = nothing, ny_core = nothing, relief::Real = 0.0)
    dxc, dyc = max(1.0, dx_core), max(1.0, dy_core)
    nx_pad, ny_pad = max(0, nx_pad), max(0, ny_pad)
    cx_mid  = (minimum(sx) + maximum(sx)) / 2
    cy_mid  = (minimum(sy) + maximum(sy)) / 2
    nx_core = nx_core === nothing ? max(1, ceil(Int, max(maximum(sx) - minimum(sx), dxc) / dxc)) : max(1, Int(nx_core))
    ny_core = ny_core === nothing ? max(1, ceil(Int, max(maximum(sy) - minimum(sy), dyc) / dyc)) : max(1, Int(ny_core))
    x0      = cx_mid - nx_core * dxc / 2
    y0      = cy_mid - ny_core * dyc / 2

    padx = [dxc * pad_factor^i for i in 1:nx_pad]
    pady = [dyc * pad_factor^i for i in 1:ny_pad]
    dx   = vcat(reverse(padx), fill(dxc, nx_core), padx)
    dy   = vcat(reverse(pady), fill(dyc, ny_core), pady)
    origin_x, origin_y = x0 - sum(padx), y0 - sum(pady)

    # uniform layers through the relief, then geometric growth to the target depth
    depth_target = depth_mult * min(nx_core * dxc, ny_core * dyc)
    dz = Float64[]
    z = 0.0
    while z < Float64(relief) - 1e-9
        push!(dz, z_first); z += z_first
    end
    t = Float64(z_first)
    while z < relief + depth_target
        push!(dz, t); z += t; t *= z_factor
    end

    x_edges = origin_x .+ vcat(0.0, cumsum(dx))
    y_edges = origin_y .+ vcat(0.0, cumsum(dy))
    z_edges = vcat(0.0, cumsum(dz))
    centres(e) = (e[1:end-1] .+ e[2:end]) ./ 2
    pad_extent = min(sum(padx), sum(pady))
    (dx = dx, dy = dy, dz = dz, origin = [origin_x, origin_y, 0.0], rho_bg = Float64(rho_bg),
     nx = length(dx), ny = length(dy), nz = length(dz),
     nx_core = nx_core, ny_core = ny_core, nx_pad = nx_pad, ny_pad = ny_pad,
     dx_core = dxc, dy_core = dyc, z_first = Float64(z_first),
     x_edges = x_edges, y_edges = y_edges, z_edges = z_edges,
     x_centers = centres(x_edges), y_centers = centres(y_edges), z_centers = centres(z_edges),
     depth = z_edges[end] - relief, pad_extent = pad_extent,
     pad_ok = pad_extent >= z_edges[end] - relief,
     depth_ok = z_edges[end] - relief >= 0.3 * min(nx_core * dxc, ny_core * dyc))
end

#----- sensors must stay above the ground -----

# deepen the ground of any column whose sensor would sit inside an earth cell, so every sensor is in air
function _apply_sensor_constraints(m, sx, sy, sdepth, ground)
    adjusted    = copy(ground)
    adjustments = NamedTuple[]
    for i in eachindex(sx)
        ix = searchsortedlast(m.x_edges, sx[i]); iy = searchsortedlast(m.y_edges, sy[i])
        (1 <= ix <= m.nx && 1 <= iy <= m.ny) || continue
        # the layer holding the sensor must be air, so the ground drops to its bottom
        k = searchsortedlast(m.z_edges, sdepth[i] - 1e-6)
        1 <= k <= m.nz || continue
        floor_depth = m.z_edges[k + 1]
        if adjusted[ix, iy] < floor_depth - 1e-9
            push!(adjustments, (datum = i, ix = ix, iy = iy, sensor_depth = sdepth[i],
                                ground_from = adjusted[ix, iy], ground_to = floor_depth))
            adjusted[ix, iy] = floor_depth
        end
    end
    adjusted, adjustments
end

# background density in the earth, NaN above the ground
function _build_density_model(m, ground)
    A = fill(m.rho_bg, m.nx, m.ny, m.nz)
    @inbounds for ix in 1:m.nx, iy in 1:m.ny, iz in 1:m.nz
        m.z_centers[iz] < ground[ix, iy] && (A[ix, iy, iz] = NaN)
    end
    A
end

"""
    build_mesh_bundle(sx, sy, elevation; rho_bg, dx_core, dy_core, nx_core, ny_core, nx_pad, ny_pad,
                      pad_factor, z_first, z_factor, depth_mult, cov_value, topography=nothing,
                      surface_elevation=0.0)

Mesh, start model and covariance mask for sensors at north `sx`, east `sy` and `elevation`.
"""
function build_mesh_bundle(sx, sy, elevation; rho_bg, dx_core, dy_core, nx_core, ny_core, nx_pad, ny_pad,
                           pad_factor, z_first, z_factor, depth_mult, cov_value,
                           topography = nothing, surface_elevation::Real = 0.0, invert_padding::Bool = true)
    args = (; rho_bg, dx_core, dy_core, nx_core, ny_core, nx_pad, ny_pad, pad_factor, z_first, z_factor, depth_mult)
    base = design_mesh(sx, sy; args..., relief = 0.0)
    elev = topography === nothing ? fill(Float64(surface_elevation), base.nx, base.ny) :
           [Float64(topography(x, y)) for x in base.x_centers, y in base.y_centers]
    all(isfinite, elev) || error("the topography returned a non-finite elevation")
    top_elevation = maximum(elev)
    if topography !== nothing
        # relief layers stay uniform from the highest ground down to the lowest
        base = design_mesh(sx, sy; args..., relief = top_elevation - minimum(elev))
    end
    ground = top_elevation .- elev                                  # depth of the ground per column
    ground, adjustments = _apply_sensor_constraints(base, sx, sy, top_elevation .- elevation, ground)
    A    = _build_density_model(base, ground)
    mask = Int.(.!isnan.(A))
    if !invert_padding
        # the lateral padding is held fixed; the layers are the model and all stay free
        keep = falses(size(mask))
        keep[(nx_pad + 1):(base.nx - nx_pad), (ny_pad + 1):(base.ny - ny_pad), :] .= true
        mask[.!keep] .= 0
    end
    (; base..., A, top_elevation, cov_mask = mask, cov_value = Float64(cov_value),
       topo_active = topography !== nothing,
       ground_min = minimum(elev), ground_max = top_elevation,
       sensor_columns = (nadjusted = length(adjustments), adjustments = adjustments),
       n_air = count(isnan, A),
       gbytes = 4.0 * length(sx) * count(!isnan, A) / 2^30)
end

function _print_sensor_report(m)
    m.sensor_columns.nadjusted == 0 && return nothing
    println(@sprintf("Ground lowered in %d column(s) so that every sensor sits in air:", m.sensor_columns.nadjusted))
    for row in m.sensor_columns.adjustments[1:min(end, 10)]
        println(@sprintf("  datum %d: xcell=%d ycell=%d sensor depth %.1f m, ground %.1f -> %.1f m",
                         row.datum, row.ix, row.iy, row.sensor_depth, row.ground_from, row.ground_to))
    end
    m.sensor_columns.nadjusted > 10 && println(@sprintf("  ... %d more column(s) omitted.", m.sensor_columns.nadjusted - 10))
    nothing
end

"""
    MakeMesh3D(data_file; out_model, out_covariance, topography=nothing, surface_elevation=nothing,
               background_density=0.0, cell_width_frac=1.0, n_pad=6, pad_factor=1.4, first_layer_div=2.0,
               vertical_factor=1.15, depth_mult=0.5, cov_value=0.3, cov_apply=2, verbose=true) -> bundle

Design an inversion mesh around a gravity survey and write the start model (`mesh_start_model.den`, usable as the prior too) and the covariance (`C3.dat`) next to the data.
"""
function MakeMesh3D(data_file::AbstractString;
                    out_model::AbstractString = joinpath(dirname(abspath(data_file)), "mesh_start_model.den"),
                    out_covariance::AbstractString = joinpath(dirname(out_model), "C3.dat"),
                    topography = nothing,
                    surface_elevation::Union{Nothing, Real} = nothing,
                    background_density::Real = 0.0,
                    cov_value::Real = 0.3,
                    cov_apply::Int = 2,
                    invert_padding::Bool = true,
                    cell_width_frac::Real = 1.0,
                    n_pad::Int = 6,
                    pad_factor::Real = 1.4,
                    first_layer_div::Real = 2.0,
                    vertical_factor::Real = 1.15,
                    depth_mult::Real = 0.5,
                    verbose::Bool = true)
    isfile(data_file) || error("data file not found: $data_file")
    d = ReadData3D(data_file)
    sx, sy, elevation = d.x, d.y, d.elevation
    surface0 = surface_elevation === nothing ? minimum(elevation) : Float64(surface_elevation)

    spacing = data_spacing(sx, sy)
    dx_core = _snap(max(spacing.mean, spacing.median) * cell_width_frac, 5:5:10000)
    nx_core = max(1, ceil(Int, max(maximum(sx) - minimum(sx), 1.0) / dx_core))
    ny_core = max(1, ceil(Int, max(maximum(sy) - minimum(sy), 1.0) / dx_core))
    z_first = _snap(dx_core / first_layer_div, 1:1:5000)

    m = build_mesh_bundle(sx, sy, elevation; rho_bg = background_density, dx_core, dy_core = dx_core,
                          nx_core, ny_core, nx_pad = n_pad, ny_pad = n_pad, pad_factor, z_first,
                          z_factor = vertical_factor, depth_mult, cov_value, topography,
                          surface_elevation = surface0, invert_padding)

    if verbose
        @printf("%d data, spacing %.0f m mean, %.0f m nearest; sensors %.0f to %.0f m; ground %s\n",
                ndata(d), spacing.mean, spacing.median, minimum(elevation), maximum(elevation),
                m.topo_active ? @sprintf("%.0f to %.0f m", m.ground_min, m.ground_max) :
                                @sprintf("flat at %.0f m", surface0))
        @printf("grid %d x %d x %d (%d cells, %d air); core %d x %d @ %.0f m; first layer %.0f m; pad %.2f km/side; depth %.2f km; G ~ %.2f GB\n",
                m.nx, m.ny, m.nz, m.nx * m.ny * m.nz, m.n_air, m.nx_core, m.ny_core, m.dx_core, m.z_first,
                m.pad_extent / 1000, m.depth / 1000, m.gbytes)
        m.pad_ok || @warn "the lateral padding is narrower than the model depth"
        m.depth_ok || @warn "the model is shallow compared with the core width"
    end
    model = GravModel3D(m.dx, m.dy, m.dz, m.A, m.origin; top_elevation = m.top_elevation)
    WriteModel3D(out_model, model; title = "MakeMesh3D start model")
    WriteCov3D(out_covariance, Cov3D(m.cov_mask; smoothing = m.cov_value, n_smooth = cov_apply))
    _print_sensor_report(m)
    verbose && (@printf("wrote %s\n", out_model); @printf("wrote %s\n", out_covariance))
    m
end
