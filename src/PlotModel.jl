# Shared plotting utilities
# Author: @pankajkmishra
# Style, colour ranges, slice extraction, body outlines and the topography staircase, reused by the data,
# slice and run plots. Core grid/index/depth utilities (edges_from_centers, core_indices,
# z_indices_for_max_depth) are in CoreUtils3D.jl and always available.
# Map axes are labelled in degrees whenever the CRS can be transformed: the data stay in projected metres and
# only the ticks are placed at round longitudes and latitudes, so cells and distances are undistorted.

using CairoMakie

#---------- style ----------

const _GRAV_PX_PER_UNIT = 3
const _GRAV_COLORMAP = Reverse(:Spectral)

# the topography, the body outlines and the station markers are all drawn as thin black lines
const _GRAV_TOPO_WIDTH     = 1.0
const _GRAV_OUTLINE_WIDTH  = 1.0
const _GRAV_STATION_STROKE = 0.4

# heavier spines and ticks, no grid, regular-weight titles
const _GRAV_AX = (spinewidth = 2.0, xtickwidth = 2.0, ytickwidth = 2.0, xticksize = 7, yticksize = 7,
                  xgridvisible = false, ygridvisible = false, xlabelfont = :regular, ylabelfont = :regular,
                  titlefont = :regular)

_grav_axis(pos; kwargs...) = Axis(pos; _GRAV_AX..., kwargs...)

_grav_save(path, figure) =
    (mkpath(dirname(abspath(path))); save(path, figure; px_per_unit = _GRAV_PX_PER_UNIT, backend = CairoMakie); String(path))

# (a), (b), ... in the top-left corner of the given layout positions
_grav_labels!(fig, positions) = for (n, p) in enumerate(positions)
    Label(fig[p..., TopLeft()], "($('a' + n - 1))"; fontsize = 22, padding = (0, 12, 10, 0), halign = :left)
end

#---------- colour range ----------

"""
    compute_colorrange(values; density_range=nothing, symmetric=true) -> (cmin, cmax)

Colour range of a model: the given range, or the 2–98 % quantiles of the finite values, made symmetric about zero for a signed density contrast.
"""
function compute_colorrange(R::AbstractArray; density_range = nothing, symmetric::Bool = true)
    if density_range === nothing
        vals = filter(isfinite, vec(R))
        isempty(vals) && (vals = [0.0, 1.0])
        cmin, cmax = quantile(vals, (0.02, 0.98))
        if symmetric
            a = max(abs(cmin), abs(cmax))
            cmin, cmax = -a, a
        end
    else
        cmin, cmax = minmax(Float64(density_range[1]), Float64(density_range[2]))
        (isfinite(cmin) && isfinite(cmax)) || error("density_range must contain finite numbers")
    end
    if cmin == cmax
        e = max(1e-12, 1e-6 * abs(cmin))
        cmin -= e; cmax += e
    end
    cmin, cmax
end

# a symmetric range from the strongest contrasts, so a body in a few cells is not flattened by the quantiles
function _body_density_range(R)
    vals = filter(isfinite, vec(R))
    a    = isempty(vals) ? 1.0 : max(quantile(abs.(vals), 0.995), 1e-6)
    (-a, a)
end

#---------- model arrays ----------

"""
    prepare_model_arrays(model; with_padding=false, max_depth=nothing, pad_tol=0.2)

Cell centres, values and index ranges of the part of a model to plot; `max_depth` (m) is measured from the top of the mesh.
"""
function prepare_model_arrays(M::GravModel3D; with_padding::Bool = false,
                              max_depth::Union{Nothing, Real} = nothing, pad_tol::Real = 0.2)
    ix, iy = with_padding ? (1:length(M.cx), 1:length(M.cy)) : lateral_core_ranges(M; tol = pad_tol)
    kz = max_depth === nothing ? (1:length(M.cz)) : z_indices_for_max_depth(M.cz, float(max_depth))
    M.cx[ix], M.cy[iy], M.cz[kz], M.A[ix, iy, kz], ix, iy, kz
end

#---------- slices ----------

"""
    model_slice(model, view, index; with_padding=false, max_depth=nothing, pad_tol=0.2)

One slice of a model as `(h_edges, v_edges, values, title)` in metres for `heatmap!`: - `:xy`: depth slice `index`, east horizontal, north vertical - `:xz`: section at east column `index`, north horizontal, depth vertical - `:yz`: section at north row `index`, east horizontal, depth vertical
"""
function model_slice(M::GravModel3D, view::Symbol, index::Integer; kwargs...)
    x, y, z, R, ix, iy, kz = prepare_model_arrays(M; kwargs...)
    xe, ye = M.x[first(ix):last(ix)+1], M.y[first(iy):last(iy)+1]
    de = M.z[first(kz):last(kz)+1]
    if view == :xy
        k = clamp(index, 1, length(z))
        return ye, xe, permutedims(R[:, :, k]), @sprintf("Depth %.0f to %.0f m", M.z[kz[k]], M.z[kz[k]+1])
    elseif view == :xz
        j = clamp(index, 1, length(y))
        return xe, de, R[:, j, :], @sprintf("East %.0f m", y[j])
    elseif view == :yz
        i = clamp(index, 1, length(x))
        return ye, de, R[i, :, :], @sprintf("North %.0f m", x[i])
    end
    throw(ArgumentError("view must be :xy, :xz or :yz"))
end

_slice_count(M, view; kwargs...) =
    (a = prepare_model_arrays(M; kwargs...); view == :xy ? length(a[3]) : view == :xz ? length(a[2]) : length(a[1]))

#---------- geographic axes ----------

# a round degree step near the asked-for spacing, so the labels land on 1, 2 or 5 of a decade
function _nice_degree_step(target::Real)
    target <= 0 && return 1.0
    p = 10.0^floor(log10(target))
    for m in (1.0, 2.0, 5.0)
        m * p >= target && return m * p
    end
    10p
end

# a degree value with its hemisphere, to as many decimals as the tick step needs
_format_degree(g::Real, dir::Symbol, digits::Integer) =
    string(round(abs(g), digits = digits), "°", dir == :east ? (g < 0 ? "W" : "E") : (g < 0 ? "S" : "N"))

# ticks at round longitudes or latitudes for an axis in projected metres, nothing when the CRS cannot be transformed
# the degrees of one axis depend slightly on the other, so both are taken through the reference point of the map
function _geo_tickfinder(crs, dir::Symbol, east_ref::Real, north_ref::Real)
    fwd = _geo_transform(crs)
    inv = _geo_transform(crs; inverse = true)
    (fwd === nothing || inv === nothing) && return nothing
    centre     = fwd((east_ref, north_ref))
    lon0, lat0 = centre[1], centre[2]
    (isfinite(lon0) && isfinite(lat0)) || return nothing
    to_degree(v) = dir == :east ? fwd((v, north_ref))[1] : fwd((east_ref, v))[2]
    to_metre(g)  = dir == :east ? inv((g, lat0))[1] : inv((lon0, g))[2]
    function (vmin, vmax)
        lo, hi = minmax(to_degree(vmin), to_degree(vmax))
        (isfinite(lo) && isfinite(hi) && hi > lo) || return (Float64[], String[])
        step    = _nice_degree_step((hi - lo) / 4)
        digits  = clamp(ceil(Int, -log10(step)), 0, 6)
        degrees = (ceil(lo / step) * step):step:hi
        isempty(degrees) && return (Float64[], String[])
        (Float64[to_metre(g) for g in degrees], [_format_degree(g, dir, digits) for g in degrees])
    end
end

# label and ticks of one map axis: degrees when the CRS can be transformed, projected metres otherwise
function _geo_axis(crs, dir::Symbol, east_ref::Real, north_ref::Real)
    ticks = _geo_tickfinder(crs, dir, east_ref, north_ref)
    ticks === nothing || return (label = dir == :east ? "Longitude" : "Latitude", ticks = ticks)
    metres = _is_local_crs(crs) ? (dir == :east ? "Y East (m)" : "X North (m)") :
                                  (dir == :east ? "Easting (m)" : "Northing (m)")
    (label = metres, ticks = nothing)
end

# the horizontal and vertical axes of a slice: :xy is east over north, :xz north over depth, :yz east over depth
function _slice_axes(view::Symbol, crs, east_ref::Real, north_ref::Real)
    east  = _geo_axis(crs, :east, east_ref, north_ref)
    north = _geo_axis(crs, :north, east_ref, north_ref)
    depth = (label = "Depth (m)", ticks = nothing)
    view == :xy ? (east, north) : view == :xz ? (north, depth) : (east, depth)
end

# axis keywords for a horizontal and a vertical axis spec, leaving the ticks automatic where there are none
_axis_kw(h, v) = merge((; xlabel = h.label, ylabel = v.label),
                       h.ticks === nothing ? NamedTuple() : (; xticks = h.ticks),
                       v.ticks === nothing ? NamedTuple() : (; yticks = v.ticks))

# the centre of a model's mesh as (east, north), the reference point its degree ticks are taken through
_model_centre(M::GravModel3D) = ((M.y[1] + M.y[end]) / 2, (M.x[1] + M.x[end]) / 2)

# the centre of a survey as (east, north)
_data_centre(d) = ((minimum(d.y) + maximum(d.y)) / 2, (minimum(d.x) + maximum(d.x)) / 2)

#---------- outlines ----------

# cell-edge outline of the cells where `mask` is true, as line segments for `linesegments!`
function _cell_outline(mask::AbstractMatrix{Bool}, he::AbstractVector, ve::AbstractVector)
    seg = Point2f[]
    nh, nv = size(mask)
    inside(i, j) = 1 <= i <= nh && 1 <= j <= nv && mask[i, j]
    for i in 1:nh, j in 1:nv
        mask[i, j] || continue
        h1, h2, v1, v2 = he[i], he[i+1], ve[j], ve[j+1]
        inside(i - 1, j) || push!(seg, Point2f(h1, v1), Point2f(h1, v2))
        inside(i + 1, j) || push!(seg, Point2f(h2, v1), Point2f(h2, v2))
        inside(i, j - 1) || push!(seg, Point2f(h1, v1), Point2f(h2, v1))
        inside(i, j + 1) || push!(seg, Point2f(h1, v2), Point2f(h2, v2))
    end
    seg
end

# the top of the first earth cell of each column of a section, as a staircase
function _topography_steps(values::AbstractMatrix, he::AbstractVector, ve::AbstractVector)
    hs, vs = Float64[], Float64[]
    for i in axes(values, 1)
        k = findfirst(!isnan, view(values, i, :))
        k === nothing && continue
        push!(hs, he[i], he[i+1]); push!(vs, ve[k], ve[k])
    end
    hs, vs
end
