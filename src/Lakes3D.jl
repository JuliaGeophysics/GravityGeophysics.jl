# Lakes for the 3D mesh: water between a lake's level and its bed, as in MTGeophysics (mask 9)
# Author: @pankajkmishra
# A lake is `(polygon = rings, level = m a.s.l.)`, optionally with `depth` (m, deepest point below the level) and
# `taper` (m): rings are closed or open lists of (x, y) vertices in `lakes_crs` (default WGS84 lon, lat; "frame" for
# the mesh frame, x north and y east in m, as MTGeophysics' 2D lakes); one ring is
# enough, several (outer shores, islands, other parts) are combined by the even-odd rule, so an island is dry land.
# MakeMesh3D takes the ground surface (DEM or `surface`) as the lake bed, as MT does; a DEM shows the water
# surface over a lake, so a lake may instead carry a modelled bed: `depth` below the level, reached `taper` metres
# from the shore (0 = a flat bottom at `depth`), never above the ground. Cells between the level and the bed are
# water: at the water contrast 1000 − ρ_ref, fixed, mask 9, and outside the smoothing. Stations may fly over a lake, not below
# its level. `read_lakes` builds the list from a shapefile.

using GeoInterface
using Shapefile
import Proj

#---------- lake list ----------

# rings as Vector{Vector{NTuple{2,Float64}}}: a single ring may be given as a vector of points
_lake_rings(p) = first(p) isa Union{Tuple, AbstractVector{<:Real}} ? [[(Float64(v[1]), Float64(v[2])) for v in p]] :
                                                                      [[(Float64(v[1]), Float64(v[2])) for v in r] for r in p]

function _check_lake(l)
    hasproperty(l, :polygon) && hasproperty(l, :level) ||
        throw(ArgumentError("a lake needs `polygon` and `level` (m above sea level)"))
    isfinite(l.level) || throw(ArgumentError("lake level must be finite"))
    depth = hasproperty(l, :depth) ? Float64(l.depth) : NaN
    taper = hasproperty(l, :taper) ? Float64(l.taper) : 0.0
    (isnan(depth) || depth > 0) && taper >= 0 || throw(ArgumentError("lake depth must be positive and taper nonnegative"))
    rings = _lake_rings(l.polygon)
    all(r -> length(r) >= 3, rings) || throw(ArgumentError("a lake ring needs at least three vertices"))
    (; rings, level = Float64(l.level), depth, taper)
end

# lake rings from `lakes_crs` into the mesh frame (x north, y east, m), through the same local TM as the DEM sampling
function _lakes_to_frame(d::Data, lakes, lakes_crs::AbstractString)
    isempty(lakes) && return NamedTuple[]
    frame = lowercase(strip(lakes_crs)) == "frame"
    ctx = frame ? nothing : _local_tm_to_target_context(d, lakes_crs)
    back = frame ? nothing : Proj.inv(ctx.trans)
    map(lakes) do l
        k = _check_lake(l)
        rings = frame ? k.rings :
            [[_target_xy_to_local_xy(ctx, back, x, y) for (x, y) in r] for r in k.rings]
        xs = [p[1] for r in rings for p in r]; ys = [p[2] for r in rings for p in r]
        (; rings, level = k.level, depth = k.depth, taper = k.taper, box = (extrema(xs), extrema(ys)))
    end
end

#---------- geometry ----------

# even-odd point in rings (closed or not)
function _in_rings(x, y, rings)
    inside = false
    for r in rings, i in eachindex(r)
        (x1, y1), (x2, y2) = r[i], r[mod1(i + 1, length(r))]
        if (y1 > y) != (y2 > y) && x < x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            inside = !inside
        end
    end
    inside
end

# distance from (x, y) to the nearest ring edge
function _ring_distance(x, y, rings)
    dmin = Inf
    for r in rings, i in eachindex(r)
        (x1, y1), (x2, y2) = r[i], r[mod1(i + 1, length(r))]
        dx, dy = x2 - x1, y2 - y1
        t = dx == dy == 0 ? 0.0 : clamp(((x - x1) * dx + (y - y1) * dy) / (dx^2 + dy^2), 0, 1)
        dmin = min(dmin, hypot(x - x1 - t * dx, y - y1 - t * dy))
    end
    dmin
end

# per column: water surface z (down, sea level 0; NaN off the lakes) and the lake bed laid into `surface` (z of the ground)
function _apply_lakes(x_centers, y_centers, surface::AbstractMatrix, lakes)
    level = fill(NaN, size(surface))
    bed = copy(surface)
    for l in lakes, (i, x) in enumerate(x_centers), (j, y) in enumerate(y_centers)
        (l.box[1][1] <= x <= l.box[1][2] && l.box[2][1] <= y <= l.box[2][2] && _in_rings(x, y, l.rings)) || continue
        level[i, j] = -l.level
        isnan(l.depth) && continue
        f = l.taper > 0 ? min(1.0, _ring_distance(x, y, l.rings) / l.taper) : 1.0
        bed[i, j] = max(bed[i, j], -l.level + l.depth * f)
    end
    level, bed
end

#---------- shapefile ----------

_points(ring) = [(Float64(GeoInterface.x(p)), Float64(GeoInterface.y(p))) for p in GeoInterface.getpoint(ring)]
function _geometry_rings(g)
    t = GeoInterface.geomtrait(g)
    t isa GeoInterface.PolygonTrait && return [_points(r) for r in GeoInterface.getring(g)]
    t isa GeoInterface.MultiPolygonTrait && return [_points(r) for p in GeoInterface.getgeom(g) for r in GeoInterface.getring(p)]
    throw(ArgumentError("lake geometries must be polygons, got $(typeof(t))"))
end

"""
    read_lakes(path; level, select=nothing, crs=nothing, depth=nothing, taper=0, topo_file="", topo_crs="EPSG:4326")

Lakes for `MakeMesh3D` from a polygon shapefile, as `(polygon, level[, depth, taper])` in WGS84 (lon, lat), the
default `lakes_crs`. Islands and multi-part lakes are kept as rings.
- `level`: water level (m a.s.l.): a number for every lake, the name of an attribute, or `:dem` for the median of
  the DEM `topo_file` (in `topo_crs`) over the lake; a DEM is flat over lakes, at the water level
- `select`: keep the rows for which `select(row)` is true (e.g. `r -> r.jarviNro == "04.411.1.001"`)
- `crs`: CRS of the shapefile; nothing reads it from the `.prj` next to it
- `depth`, `taper`: modelled bed (see `Lakes3D.jl`): a number, an attribute name or a function of the row; nothing
  keeps the ground as the bed
"""
function read_lakes(path::AbstractString; level, select = nothing, crs = nothing, depth = nothing, taper = 0.0,
                    topo_file::AbstractString = "", topo_crs::AbstractString = "EPSG:4326")
    isfile(path) || error("lake shapefile not found: $path")
    src = crs !== nothing ? crs : (wkt = last(detect_shapefile_crs(path)); isempty(wkt) ? "EPSG:4326" : wkt)
    tr = Proj.Transformation(src, "EPSG:4326"; always_xy = true)
    level === :dem && (isfile(topo_file) || error("level = :dem needs the DEM as topo_file"))
    dem = level === :dem ? _ensure_topography(topo_file) : nothing
    to_dem = level === :dem ? Proj.Transformation(src, topo_crs; always_xy = true) : nothing
    value(spec, row) = spec === nothing ? nothing : spec isa Real ? Float64(spec) : spec isa Function ? Float64(spec(row)) :
                                         Float64(getproperty(row, Symbol(spec)))
    lakes = NamedTuple[]
    for row in Shapefile.Table(path)
        select === nothing || select(row) || continue
        g = GeoInterface.geometry(row)
        g === missing && continue
        rings = _geometry_rings(g)
        h = level === :dem ? _dem_level(dem, to_dem, rings) : value(level, row)
        isfinite(h) || error("no water level for a lake in $path")
        lake = (; polygon = [[(p = tr(v); (Float64(p[1]), Float64(p[2]))) for v in r] for r in rings], level = h)
        z = value(depth, row)
        push!(lakes, z === nothing ? lake : (; lake..., depth = z, taper = value(taper, row)))
    end
    isempty(lakes) && error("no lake selected in $path")
    lakes
end

# median of the DEM on a grid of up to 40 × 40 points inside the rings (the vertices for a lake smaller than that)
function _dem_level(dem, to_dem, rings)
    xs = [p[1] for r in rings for p in r]; ys = [p[2] for r in rings for p in r]
    (x0, x1), (y0, y1) = extrema(xs), extrema(ys)
    pts = [(x, y) for x in range(x0, x1; length = 40), y in range(y0, y1; length = 40) if _in_rings(x, y, rings)]
    isempty(pts) && (pts = [p for r in rings for p in r])
    q = [to_dem(p) for p in pts]
    all(p -> _in_dem(dem, p[1], p[2]), q) || error("a lake lies outside the DEM; its level cannot be read from it")
    v = filter(isfinite, [_sample_topography(dem, p[1], p[2]) for p in q])
    isempty(v) ? NaN : median(v)
end
