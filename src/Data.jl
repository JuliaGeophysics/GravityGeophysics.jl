# 3D gravity data files
# Author: @pankajkmishra
# Reads and writes the .xyz survey layout, extended with an error column and the geographic position
# Internally x north, y east and the sensor elevation; depths are measured down from the top of the model mesh,
# the highest topography, so a sensor's depth is `top_elevation - elevation` (negative above the mesh top)
# Every station also carries its latitude and longitude, derived from the projected position and the CRS header
# when the file does not give them, so a file is georeferenced on its own and the plots can label in degrees

using Printf
using Proj

#---------- data ----------

"""
    GravData3D

A gravity survey, one row per datum.
"""
Base.@kwdef mutable struct GravData3D
    x::Vector{Float64}
    y::Vector{Float64}
    elevation::Vector{Float64}
    value::Vector{Float64}
    error::Vector{Float64}
    lat::Vector{Float64} = Float64[]
    lon::Vector{Float64} = Float64[]
    crs::String = "local"
    title::String = "Vertical gravity anomaly"
end

"""
    ndata(data) -> Int

Number of stations in a survey.
"""
ndata(d::GravData3D) = length(d.value)

function _validate_data(d::GravData3D)
    n = length(d.value)
    length(d.x) == length(d.y) == length(d.elevation) == length(d.error) == n ||
        throw(DimensionMismatch("coordinate, value and error columns must have the same length"))
    n > 0 || throw(ArgumentError("no data"))
    all(isfinite, d.x) && all(isfinite, d.y) && all(isfinite, d.elevation) || throw(ArgumentError("coordinates must be finite"))
    all(isfinite, d.value) || throw(ArgumentError("data values must be finite"))
    isempty(d.lat) || length(d.lat) == n || throw(DimensionMismatch("the latitude column must cover every station"))
    isempty(d.lon) || length(d.lon) == n || throw(DimensionMismatch("the longitude column must cover every station"))
    d
end

"""
    has_geographic(data) -> Bool

True when a survey carries a latitude and a longitude for every station.
"""
has_geographic(d::GravData3D) = length(d.lat) == length(d.lon) == ndata(d) && ndata(d) > 0

#---------- header keys ----------

# "# Key : value" lines of the data header, matched like the control files
const _DATA_HEADER_SPEC = (
    "CRS" => :crs,
)

_header_key(s) = lowercase(join(split(strip(s)), " "))

#---------- CRS ----------

_is_local_crs(crs) = lowercase(strip(crs)) in ("", "local", "none", "model")

function _is_geographic_crs(crs)
    _is_local_crs(crs) && return false
    c = uppercase(strip(crs))
    c in ("EPSG:4326", "WGS84", "EPSG:4258", "EPSG:4269", "CRS84", "OGC:CRS84")
end

# WGS84 UTM zone of a longitude/latitude centre
function _utm_crs(lon0::Real, lat0::Real)
    zone = clamp(floor(Int, (lon0 + 180) / 6) + 1, 1, 60)
    lat0 >= 0 ? "EPSG:$(32600 + zone)" : "EPSG:$(32700 + zone)"
end

const _WGS84 = "EPSG:4326"

# a transformation between a projected CRS and longitude/latitude, or nothing when the CRS is local or unknown
function _geo_transform(crs; inverse::Bool = false)
    _is_local_crs(crs) && return nothing
    try
        inverse ? Proj.Transformation(_WGS84, String(crs); always_xy = true) :
                  Proj.Transformation(String(crs), _WGS84; always_xy = true)
    catch err
        @warn "'$crs' is not a CRS this build of Proj can transform, positions stay in metres" err
        nothing
    end
end

# latitude and longitude of projected eastings and northings, empty vectors when the CRS cannot be transformed
function _geographic_columns(crs, east::AbstractVector, north::AbstractVector)
    trans = _geo_transform(crs)
    trans === nothing && return (Float64[], Float64[])
    lat, lon = similar(Vector{Float64}, length(east)), similar(Vector{Float64}, length(east))
    for i in eachindex(east)
        p = trans((east[i], north[i]))
        lon[i], lat[i] = p[1], p[2]
    end
    (lat, lon)
end

#---------- reading ----------

"""
    ReadData3D(path; crs=nothing) -> GravData3D

Read a gravity survey in the `.xyz` layout: `easting northing elevation value` records (gz, mGal), optionally followed by `error` (mGal) and by `latitude longitude` (degrees), after `#` header lines.
"""
function ReadData3D(path::AbstractString; crs::Union{Nothing, AbstractString} = nothing)
    isfile(path) || error("data file not found: $path")
    lookup = Dict(_header_key(k) => v for (k, v) in _DATA_HEADER_SPEC)
    head   = Dict{Symbol, String}()
    title  = ""
    cols   = [Float64[] for _ in 1:7]
    ncol   = 0
    for (n, raw) in enumerate(eachline(path))
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, '#')
            body = strip(line[2:end])
            isempty(title) && !occursin(':', body) && (title = body)
            occursin(':', body) || continue
            key, value = strip.(split(body, ':'; limit = 2))
            field = get(lookup, _header_key(key), nothing)
            field === nothing && continue
            haskey(head, field) && error("$path:$n: duplicate header key '$key'")
            head[field] = value
            continue
        end
        tokens = split(first(split(line, '#')))
        length(tokens) in (4, 5, 7) ||
            error("$path:$n: expected 4, 5 or 7 columns, got $(length(tokens))")
        ncol == 0 && (ncol = length(tokens))
        length(tokens) == ncol || error("$path:$n: $(length(tokens)) columns where earlier records have $ncol")
        for k in 1:ncol
            v = tryparse(Float64, replace(tokens[k], 'D' => 'E', 'd' => 'e'))
            v === nothing && error("$path:$n: non-numeric value '$(tokens[k])'")
            push!(cols[k], v)
        end
    end
    ncol == 0 && error("$path: no data records")
    errors = if ncol >= 5
        cols[5]
    else
        @warn "$path has no error column; every error is 1 mGal, set 'Error floor (mGal)' in the inversion control"
        ones(length(cols[1]))
    end
    all(e -> isfinite(e) && e > 0, errors) || error("$path: errors must be positive and finite")

    source      = strip(crs === nothing ? get(head, :crs, "local") : String(crs))
    east, north = cols[1], cols[2]
    lat, lon    = ncol == 7 ? (cols[6], cols[7]) : (Float64[], Float64[])
    if _is_geographic_crs(source)
        isempty(lat) && ((lat, lon) = (copy(north), copy(east)))
        target = _utm_crs(sum(east) / length(east), sum(north) / length(north))
        trans  = Proj.Transformation(source, target; always_xy = true)
        for i in eachindex(east)
            p = trans((east[i], north[i]))
            east[i], north[i] = p[1], p[2]
        end
        @info "Projected $(basename(path)) from $source to $target"
        source = target
    end
    # a file that does not carry the geographic columns gets them from its CRS, so every survey is georeferenced
    isempty(lat) && ((lat, lon) = _geographic_columns(source, east, north))

    _validate_data(GravData3D(
        x = north, y = east, elevation = cols[3], value = cols[4], error = errors, lat = lat, lon = lon,
        crs = String(source), title = isempty(title) ? "Gravity data" : title))
end

#---------- writing ----------

"""
    WriteData3D(path, data::GravData3D; title=data.title, write_error=true, write_geographic=true) -> path

Write a survey in the `.xyz` layout with the header keys of `ReadData3D`, the projected position first and the latitude and longitude last; a reader that wants only positions and values takes the first four columns and skips the header.
"""
function WriteData3D(path::AbstractString, d::GravData3D; title::AbstractString = d.title,
                     write_error::Bool = true, write_geographic::Bool = true)
    _validate_data(d)
    mkpath(dirname(abspath(path)))
    unit = "mGal"
    # the geographic columns ride behind the error column, so they need it even when every error is the default
    geo  = write_geographic && has_geographic(d)
    open(path, "w") do io
        println(io, "# ", title)
        println(io, "# CRS : ", d.crs)
        geo && println(io, "# Geographic CRS : ", _WGS84)
        println(io, "# easting(m) northing(m) elevation(m) value($unit)",
                write_error || geo ? " error($unit)" : "", geo ? " latitude(deg) longitude(deg)" : "")
        for i in eachindex(d.value)
            @printf(io, "%14.3f %14.3f %10.3f %14.6g", d.y[i], d.x[i], d.elevation[i], d.value[i])
            (write_error || geo) && @printf(io, " %12.6g", d.error[i])
            geo && @printf(io, " %14.8f %14.8f", d.lat[i], d.lon[i])
            println(io)
        end
    end
    String(path)
end

"""
    data_like(data, values; error=data.error, title=data.title) -> GravData3D

Copy of a survey with new values, e.g.
"""
data_like(d::GravData3D, values::AbstractVector{<:Real}; error::AbstractVector{<:Real} = d.error,
          title::AbstractString = d.title) =
    GravData3D(x = copy(d.x), y = copy(d.y), elevation = copy(d.elevation), value = Float64.(values),
               error = Float64.(error), lat = copy(d.lat), lon = copy(d.lon), crs = d.crs, title = String(title))
