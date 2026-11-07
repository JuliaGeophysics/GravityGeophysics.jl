# Gravity data: the survey frame, the native data file, CSV and Geosoft XYZ import, draping and frame checks.
# Author: @pankajkmishra
# Stations are kept in the mesh frame (x north, y east, z DOWN, metres). One data type per file.
# The native file follows the MT data file: comment lines, five '>' header lines (type, units,
# frame azimuth, origin lat lon, station count), then one row per station
#   Code  GG_Lat  GG_Lon  X(m)  Y(m)  Z(m)  Component  Value  Error
# The data frame is the mesh frame turned by `azimuth` (see rotate_data); values never change on
# a rotation. The magnetic file has one more header line (the inducing field), gravity has none.
# The survey frame takes geographic positions to mesh x (north) and y (east), in metres: without an MT
# model, a local transverse Mercator at the survey centre (the projection the resistivity plots use);
# on an existing MT model, the quadratic map from (lat, lon) to the x/y its data file holds, fitted to
# the MT sites, whatever tool made those x/y (real files are linear in lat/lon, or close to it).

#---------- survey frame ----------

struct Frame
    lat0::Float64
    lon0::Float64
    C::Union{Nothing,Matrix{Float64}}   # 6 × 2 quadratic coefficients; nothing: local TM
    rms::Float64                        # fit residual on the MT sites (m)
    latlim::NTuple{2,Float64}           # extent of the MT sites
    lonlim::NTuple{2,Float64}
    source::String
end

"""
    Frame(lat0, lon0)

Local transverse Mercator frame centred at (`lat0`, `lon0`): x north, y east, in metres.
"""
Frame(lat0::Real, lon0::Real) = Frame(lat0, lon0, nothing, 0.0, (-90.0, 90.0), (-180.0, 180.0),
                                      @sprintf("local TM at %.6f %.6f", lat0, lon0))

_tm(lat0, lon0) = Proj.Transformation("EPSG:4326",
    "+proj=tmerc +lat_0=$(Float64(lat0)) +lon_0=$(Float64(lon0)) +k=0.9996 +x_0=500000 +y_0=0 +datum=WGS84 +units=m +no_defs";
    always_xy = true)

function _tm_xy(lat0, lon0, lat::AbstractVector, lon::AbstractVector)
    tr = _tm(lat0, lon0)
    p = [tr((Float64(lo), Float64(la))) for (la, lo) in zip(lat, lon)]
    [q[2] for q in p], [q[1] - 500000.0 for q in p]
end

"""
    to_frame(f, lat, lon) -> (x, y)

Frame coordinates (m, x north, y east) of WGS84 positions.
"""
function to_frame(f::Frame, lat::AbstractVector, lon::AbstractVector)
    f.C === nothing && return _tm_xy(f.lat0, f.lon0, lat, lon)
    xy = _terms(lat .- f.lat0, lon .- f.lon0) * f.C
    xy[:, 1], xy[:, 2]
end

_terms(a, b) = hcat(ones(length(a)), a, b, a .^ 2, a .* b, b .^ 2)

"""
    north_azimuth(f, lat, lon) -> degrees

Azimuth of true north in the frame at (`lat`, `lon`), clockwise from frame x. Zero at the centre of
a plain local TM frame; an azimuth α from true north is α + north_azimuth in the frame.
"""
function north_azimuth(f::Frame, lat::Real, lon::Real)
    δ = 1e-4
    x, y = to_frame(f, [lat - δ, lat + δ], [lon, lon])
    atand(y[2] - y[1], x[2] - x[1])
end

# sites of an MT data file: origin from the fifth '>' header line, then per site code lat lon x y
function _mt_sites(path::AbstractString)
    origin = nothing; nhdr = 0
    code = String[]; lat = Float64[]; lon = Float64[]; x = Float64[]; y = Float64[]
    for l in eachline(path)
        t = strip(l)
        (isempty(t) || startswith(t, "#")) && continue
        if startswith(t, ">")
            nhdr += 1
            if nhdr == 5 && origin === nothing
                p = split(strip(t[2:end]))
                origin = parse.(Float64, p[1:2])
            end
        else
            p = split(t)
            length(p) >= 11 || continue
            p[2] in code && continue
            push!(code, p[2]); push!(lat, parse(Float64, p[3])); push!(lon, parse(Float64, p[4]))
            push!(x, parse(Float64, p[5])); push!(y, parse(Float64, p[6]))
        end
    end
    origin === nothing && throw(ArgumentError("no origin line in MT data file $path"))
    (; origin, code, lat, lon, x, y)
end

"""
    load_frame(mt_data_path) -> Frame

Frame of an existing MT model: the least-squares quadratic map from (lat, lon) to the sites' x/y
in the MT data file. The residual is reported; stations should lie within the sites' extent.
"""
function load_frame(path::AbstractString)
    s = _mt_sites(path)
    length(s.x) >= 12 || throw(ArgumentError("need at least 12 sites in $path to fit the frame"))
    M = _terms(s.lat .- s.origin[1], s.lon .- s.origin[2])
    C = M \ hcat(s.x, s.y)
    r = hcat(s.x, s.y) .- M * C
    rms = sqrt(sum(abs2, r) / length(s.x))
    @info @sprintf("frame from %d MT sites: residual rms %.2f m, max %.2f m",
                   length(s.x), rms, maximum(hypot.(r[:, 1], r[:, 2])))
    Frame(s.origin[1], s.origin[2], C, rms, extrema(s.lat), extrema(s.lon),
          @sprintf("MT sites of %s (rms %.2f m)", basename(path), rms))
end

"""
    Frame(crs, latlim, lonlim) -> Frame

Frame whose x and y are the northing and easting (m) of the projected `crs` (e.g. "EPSG:3067",
ETRS-TM35FIN), so models and maps are in that grid: the quadratic map from (lat, lon) fitted to
the CRS on a 21 × 21 grid over `latlim` × `lonlim`. The fit residual is reported; frame north is
grid north, and `north_azimuth` gives the grid convergence.
"""
function Frame(crs::AbstractString, latlim::NTuple{2,Real}, lonlim::NTuple{2,Real})
    tr = Proj.Transformation("EPSG:4326", crs; always_xy = true)
    la, lo = range(latlim...; length = 21), range(lonlim...; length = 21)
    lat, lon = [a for a in la for _ in lo], [b for _ in la for b in lo]
    p = [tr((b, a)) for (a, b) in zip(lat, lon)]
    lat0, lon0 = sum(latlim) / 2, sum(lonlim) / 2
    M = _terms(lat .- lat0, lon .- lon0)
    xy = hcat([q[2] for q in p], [q[1] for q in p])
    C = M \ xy
    rms = sqrt(sum(abs2, xy .- M * C) / length(lat))
    Frame(lat0, lon0, C, rms, Float64.(latlim), Float64.(lonlim), @sprintf("%s (quadratic fit, rms %.3f m)", crs, rms))
end

#---------- data ----------

# gz is the vertical attraction, positive down (towards a mass below the station)
const COMPONENTS = (:gz,)
const _COMPONENT_NAME = Dict(:gz => "gz")
const _COMPONENT_SYMBOL = Dict(lowercase(v) => k for (k, v) in _COMPONENT_NAME)
units(c::Symbol) = "mGal"

mutable struct Data
    x::Vector{Float64}          # mesh frame (m), z down
    y::Vector{Float64}
    z::Vector{Float64}
    lat::Vector{Float64}        # WGS84, NaN when stations were given in frame coordinates
    lon::Vector{Float64}
    value::Vector{Float64}
    error::Vector{Float64}
    code::Vector{String}        # "line:station"
    component::Symbol
    origin::Vector{Float64}     # lat, lon of the frame centre
    north::Float64              # azimuth of true north in the unrotated frame (degrees)
    azimuth::Float64            # data frame x axis, clockwise from the unrotated frame x
    rotations::Vector{Tuple{Symbol,Float64}}
    frame::String
    frame_rms::Float64          # position error of the frame fit (m)
    name::String
end

ndata(d::Data) = length(d.value)

function _check(d::Data)
    n = length(d.value)
    all(length(v) == n for v in (d.x, d.y, d.z, d.lat, d.lon, d.error, d.code)) ||
        throw(ArgumentError("data columns differ in length"))
    d.component in COMPONENTS || throw(ArgumentError("unknown component $(d.component); use one of $COMPONENTS"))
    all(isfinite, d.x) && all(isfinite, d.y) && all(isfinite, d.z) || throw(ArgumentError("nonfinite station position"))
    all(isfinite, d.value) || throw(ArgumentError("nonfinite data value"))
    all(e -> isfinite(e) && e > 0, d.error) || throw(ArgumentError("data errors must be positive and finite"))
    d
end

#---------- native data file ----------

_rotation_text(r) = join([@sprintf("%s %+.6f", k, a) for (k, a) in r], "; ")
function _parse_rotations(s::AbstractString)
    out = Tuple{Symbol,Float64}[]
    for part in split(s, ';')
        p = split(strip(part)); length(p) == 2 || continue
        push!(out, (Symbol(p[1]), parse(Float64, p[2])))
    end
    out
end

"""
    write_data_grav(path, d; value=d.value, error=d.error, comment="")

Write `d` (or predicted `value`s on its stations) as a native data file.
"""
function write_data_grav(path::AbstractString, d::Data; value::AbstractVector = d.value,
                    error::AbstractVector = d.error, comment::AbstractString = "")
    length(value) == ndata(d) == length(error) || throw(ArgumentError("value/error length differs from the data"))
    mkpath(dirname(abspath(path)))
    name = _COMPONENT_NAME[d.component]
    open(path, "w") do io
        println(io, "# GravityGeophysics data", isempty(comment) ? "" : " -- " * comment)
        @printf(io, "# north azimuth: %.6f\n", d.north)
        println(io, "# frame: ", d.frame)
        @printf(io, "# frame rms (m): %.3f\n", d.frame_rms)
        isempty(d.rotations) || println(io, "# rotated: ", _rotation_text(d.rotations))
        println(io, "# Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error")
        println(io, "> ", name)
        println(io, "> [", units(d.component), "]")
        @printf(io, "> %.6f\n", d.azimuth)
        @printf(io, "> %.6f %.6f\n", d.origin[1], d.origin[2])
        println(io, "> ", ndata(d))
        for n in 1:ndata(d)
            @printf(io, "%-16s %12.6f %12.6f %14.3f %14.3f %12.3f %-10s %15.6E %15.6E\n",
                    d.code[n], d.lat[n], d.lon[n], d.x[n], d.y[n], d.z[n], name, value[n], error[n])
        end
    end
    String(path)
end

function _load_native(path::AbstractString)
    meta = Dict{String,String}(); hdr = String[]
    rows = Vector{SubString{String}}[]
    for l in eachline(path)
        t = strip(l)
        isempty(t) && continue
        if startswith(t, "#")
            k = findfirst(':', t)
            k === nothing || (meta[lowercase(strip(t[2:k-1]))] = strip(t[k+1:end]))
        elseif startswith(t, ">")
            isempty(rows) || throw(ArgumentError("one data type per file: second header block in $path"))
            push!(hdr, strip(t[2:end]))
        else
            push!(rows, split(t))
        end
    end
    length(hdr) == 5 || throw(ArgumentError("expected five '>' header lines in $path, found $(length(hdr))"))
    component = get(_COMPONENT_SYMBOL, lowercase(hdr[1]), nothing)
    component === nothing && throw(ArgumentError("unknown data type '$(hdr[1])' in $path"))
    azimuth = parse(Float64, split(hdr[3])[1])
    origin = parse.(Float64, split(hdr[4])[1:2])
    nd = parse(Int, split(hdr[5])[1])
    nd == length(rows) || throw(ArgumentError("header says $nd stations, file has $(length(rows))"))
    all(length(r) == 9 for r in rows) || throw(ArgumentError("data rows need 9 columns in $path"))
    all(lowercase(r[7]) == lowercase(hdr[1]) for r in rows) || throw(ArgumentError("component column differs from the header in $path"))
    col(k) = [parse(Float64, r[k]) for r in rows]
    Data(col(4), col(5), col(6), col(2), col(3), col(8), col(9), [String(r[1]) for r in rows], component,
         origin, parse(Float64, get(meta, "north azimuth", "0")), azimuth,
         _parse_rotations(get(meta, "rotated", "")), get(meta, "frame", ""),
         parse(Float64, get(meta, "frame rms (m)", "0")), String(path)) |> _check
end

#---------- CSV and Geosoft XYZ import ----------

_cell(s) = String(strip(strip(s), ['"', '\'']))
_number(s) = (t = strip(s); t in ("", "*", "NaN", "nan", "NA", "-") ? NaN : parse(Float64, t))

function _read_csv(path)
    rows = [l for l in readlines(path) if !isempty(strip(l)) && !startswith(lstrip(l), "#")]
    isempty(rows) && throw(ArgumentError("empty data file $path"))
    sep = occursin(',', rows[1]) ? ',' : occursin(';', rows[1]) ? ';' : nothing
    sp(l) = _cell.(sep === nothing ? split(strip(l)) : split(l, sep))
    names = sp(rows[1])
    cells = sp.(rows[2:end])
    bad = findfirst(r -> length(r) != length(names), cells)
    bad === nothing || throw(ArgumentError("row $(bad + 1) of $path has $(length(cells[bad])) fields, header has $(length(names))"))
    names, cells, fill("", length(cells))
end

# '/' lines are comments, the last one before the data that matches the row width names the
# columns; "Line ..." and "Tie ..." lines start a line, '*' is a dummy value
function _read_xyz(path)
    heads = Vector{String}[]; cells = Vector{String}[]; lines = String[]; current = ""
    for l in eachline(path)
        t = strip(l)
        isempty(t) && continue
        if startswith(t, "/")
            isempty(cells) && push!(heads, _cell.(split(strip(t[2:end]))))
        elseif (m = match(r"^(line|tie)\s+(.*)$"i, t)) !== nothing
            current = (lowercase(m[1]) == "tie" ? "T" : "L") * replace(strip(m[2]), r"\s+" => "_")
        else
            push!(cells, _cell.(split(t))); push!(lines, current)
        end
    end
    isempty(cells) && throw(ArgumentError("no data rows in $path"))
    w = length(cells[1])
    bad = findfirst(r -> length(r) != w, cells)
    bad === nothing || throw(ArgumentError("data rows of $path differ in width"))
    k = findlast(h -> length(h) == w, heads)
    k === nothing ? String[] : heads[k], cells, lines
end

function _column(names, cells, spec, path)
    k = spec isa Integer ? Int(spec) : findfirst(n -> lowercase(n) == lowercase(String(spec)), names)
    k === nothing && throw(ArgumentError("no column '$spec' in $path; columns are $(join(names, ", "))"))
    1 <= k <= length(cells[1]) || throw(ArgumentError("column $k out of range in $path"))
    [r[k] for r in cells]
end

"""
    load_data_grav(path; component=:gz, columns=(...), crs="EPSG:4326", frame=nothing, format=:auto) -> Data

Read a native data file (the keywords are then ignored), or import CSV / Geosoft XYZ survey data.
`columns` maps roles to column names (or numbers):
- position: `lat`, `lon` (WGS84); or `easting`, `northing` in `crs`; or `north`, `east` in the frame
- `elevation` (m above sea level)
- `value`, `error` (mGal), optional `line`
Positions go to `frame` (default: local TM at the station centroid; `load_frame` for an MT model).
The values are the anomaly the density-contrast model explains, normally a complete Bouguer anomaly
at the reference density (the air cells of the mesh carry no mass). For a rotated mesh use
`rotate_data(d, m.rotation)`; `drape!` sets heights above ground; `check_frame` checks the result.
"""
function load_data_grav(path::AbstractString; component::Symbol = :gz,
                   columns::NamedTuple = NamedTuple(), crs::AbstractString = "EPSG:4326",
                   frame::Union{Nothing,Frame} = nothing, format::Symbol = :auto)
    fmt = format == :auto ? _guess_format(path) : format
    fmt == :native && return _load_native(path)
    fmt in (:csv, :xyz) || throw(ArgumentError("unknown data format $fmt"))
    names, cells, lines = fmt == :csv ? _read_csv(path) : _read_xyz(path)
    _import(path, names, cells, lines, component, columns, crs, frame)
end

function _guess_format(path)
    for l in eachline(path)
        t = strip(l)
        (isempty(t) || startswith(t, "#")) && continue
        startswith(t, ">") && return :native
        break
    end
    lowercase(splitext(path)[2]) == ".csv" ? :csv : :xyz
end

function _import(path, names, cells, lines, component, columns, crs, frame)
    component in COMPONENTS || throw(ArgumentError("unknown component $component; use one of $COMPONENTS"))
    has(k) = haskey(columns, k)
    get_num(k) = _number.(_column(names, cells, columns[k], path))
    has(:value) && has(:error) || throw(ArgumentError("columns must name `value` and `error`"))
    has(:elevation) || throw(ArgumentError("columns must name `elevation`"))
    value, err, h = get_num(:value), get_num(:error), get_num(:elevation)
    n = length(value)
    lat = fill(NaN, n); lon = fill(NaN, n)
    if has(:lat) && has(:lon)
        lat, lon = get_num(:lat), get_num(:lon)
    elseif has(:easting) && has(:northing)
        E, N = get_num(:easting), get_num(:northing)
        tr = Proj.Transformation(crs, "EPSG:4326"; always_xy = true)
        ok = isfinite.(E) .& isfinite.(N)
        for k in findall(ok)
            p = tr((E[k], N[k])); lon[k], lat[k] = p[1], p[2]
        end
    end
    framed = has(:north) && has(:east)
    framed || has(:lat) || has(:easting) ||
        throw(ArgumentError("columns must name lat/lon, easting/northing or north/east"))
    xn, ye = framed ? (get_num(:north), get_num(:east)) : (fill(0.0, n), fill(0.0, n))
    pos = framed ? (xn, ye) : (lat, lon)
    keep = isfinite.(value) .& isfinite.(err) .& isfinite.(h) .& isfinite.(pos[1]) .& isfinite.(pos[2])
    count(!, keep) > 0 && @info "$(count(!, keep)) of $n rows with dummy or missing values skipped"
    any(keep) || throw(ArgumentError("no complete rows in $path"))
    line = has(:line) ? _column(names, cells, columns[:line], path) : lines
    code = [isempty(line[k]) ? string(k) : string(line[k], ":", k) for k in 1:n]
    lat, lon, value, err, h, code, xn, ye = (v[keep] for v in (lat, lon, value, err, h, code, xn, ye))
    if framed
        x, y = xn, ye
        origin = [NaN, NaN]; north = 0.0; fdesc = "frame coordinates"; frms = 0.0
    else
        f = something(frame, Frame(sum(lat) / length(lat), sum(lon) / length(lon)))
        pad(l) = 0.1 * (l[2] - l[1])
        out = count(k -> !(f.latlim[1] - pad(f.latlim) <= lat[k] <= f.latlim[2] + pad(f.latlim) &&
                           f.lonlim[1] - pad(f.lonlim) <= lon[k] <= f.lonlim[2] + pad(f.lonlim)), eachindex(lat))
        out > 0 && @warn "$out stations lie outside the MT sites the frame was fitted to"
        x, y = to_frame(f, lat, lon)
        origin = [f.lat0, f.lon0]; fdesc = f.source; frms = f.rms
        north = north_azimuth(f, sum(lat) / length(lat), sum(lon) / length(lon))
    end
    Data(x, y, -h, lat, lon, value, err, code, component, origin, north, 0.0,
         Tuple{Symbol,Float64}[], fdesc, frms, String(path)) |> _check
end

#---------- mesh ----------

# station z on the mesh of `m`: data files keep z = −elevation, the mesh has z = 0 at m.datum (m a.s.l.)
_mesh_z(d::Data, m::Model) = d.z .+ m.datum

# column (i, j) of the mesh holding station (x, y); nothing outside the mesh
function _column_index(m::Model, x::Real, y::Real)
    (m.x[1] <= x <= m.x[end] && m.y[1] <= y <= m.y[end]) || return nothing
    (clamp(searchsortedlast(m.x, x), 1, m.nx), clamp(searchsortedlast(m.y, y), 1, m.ny))
end

"""
    drape!(d, m, height)

Put stations `height` metres (≥ 0, scalar or per station; 0 for a ground survey) above the top of
the earth cells of `m`.
"""
function drape!(d::Data, m::Model, height)
    h = height isa Real ? fill(Float64(height), ndata(d)) : Float64.(height)
    length(h) == ndata(d) || throw(ArgumentError("one height per station"))
    all(>=(0), h) || throw(ArgumentError("heights above ground must be nonnegative"))
    for n in 1:ndata(d)
        c = _column_index(m, d.x[n], d.y[n])
        c === nothing && throw(ArgumentError("station $(d.code[n]) lies outside the mesh"))
        d.z[n] = m.Z[c...] - h[n] - m.datum
    end
    d
end

"""
    check_frame(d, m)

The data frame must be the mesh frame (data azimuth = mesh rotation), and every station must lie
inside the mesh, on or above its earth cells (gz is finite and continuous on the ground, so ground
stations need no clearance). Warns when the frame fit error exceeds half a core cell.
"""
function check_frame(d::Data, m::Model)
    isapprox(rem(d.azimuth - m.rotation, 360, RoundNearest), 0; atol = 1e-9) ||
        throw(ArgumentError(@sprintf("data frame azimuth %.4f° differs from the mesh rotation %.4f°; use rotate_data(d, %.4f)",
                                     d.azimuth, m.rotation, m.rotation - d.azimuth)))
    outside = 0; below = 0
    z = _mesh_z(d, m)
    for n in 1:ndata(d)
        c = _column_index(m, d.x[n], d.y[n])
        c === nothing && (outside += 1; continue)
        z[n] <= m.Z[c...] || (below += 1)
    end
    outside == 0 || throw(ArgumentError("$outside stations lie outside the mesh"))
    below == 0 || throw(ArgumentError("$below stations lie below the top of the earth cells; drape them with drape!(d, m, height)"))
    core = min(minimum(m.dx), minimum(m.dy))
    d.frame_rms > core / 2 && @warn @sprintf("frame fit error %.1f m exceeds half the smallest cell (%.1f m)", d.frame_rms, core / 2)
    nothing
end
