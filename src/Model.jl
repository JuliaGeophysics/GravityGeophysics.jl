# Density model reader/writer
# Author: @pankajkmishra
# The block layout: nx ny nz 0 LINEAR, cell sizes, layers top down with each
# layer written row by row (y) and north to south (x reversed), then the origin and rotation
# Values are density contrasts (kg/m³); air above the topography is NaN in memory and 1e17 on file
# Depth is zero at the top of the mesh, the highest topography, and positive down; the elevation of that datum is
# written as the header line "# Top elevation (m) : value"

using Printf

#---------- model ----------

"""
    GravModel3D

A density-contrast model on a tensor mesh, x north, y east, z depth positive down.
"""
mutable struct GravModel3D
    A::Array{Float64, 3}
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    origin::Vector{Float64}
    top_elevation::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    cx::Vector{Float64}
    cy::Vector{Float64}
    cz::Vector{Float64}
    npad::NTuple{2, Int}
    name::String
end

const GRAV_AIR_VALUE = 1e17          # the air sentinel on file
const GRAV_AIR_THRESHOLD = 1e15      # values above this read as air

"""
    GravModel3D(dx, dy, dz, A, origin; top_elevation=0.0, name="")

Model from cell sizes, values and the south-west top corner `(north, east, depth)`; edges, centres and padding are derived.
"""
function GravModel3D(dx::AbstractVector{<:Real}, dy::AbstractVector{<:Real}, dz::AbstractVector{<:Real},
                     A::AbstractArray{<:Real, 3}, origin::AbstractVector{<:Real};
                     top_elevation::Real = 0.0, name::AbstractString = "")
    size(A) == (length(dx), length(dy), length(dz)) ||
        throw(DimensionMismatch("model is $(size(A)) cells, the cell sizes give $((length(dx), length(dy), length(dz)))"))
    all(>(0), dx) && all(>(0), dy) && all(>(0), dz) || throw(ArgumentError("cell sizes must be positive"))
    length(origin) == 3 || throw(ArgumentError("origin must hold north, east and depth"))
    isfinite(top_elevation) || throw(ArgumentError("top elevation must be finite"))
    x = origin[1] .+ vcat(0.0, cumsum(Float64.(dx)))
    y = origin[2] .+ vcat(0.0, cumsum(Float64.(dy)))
    z = origin[3] .+ vcat(0.0, cumsum(Float64.(dz)))
    centres(e) = (e[1:end-1] .+ e[2:end]) ./ 2
    GravModel3D(Array{Float64, 3}(A), Float64.(dx), Float64.(dy), Float64.(dz), Float64.(origin), Float64(top_elevation),
                x, y, z, centres(x), centres(y), centres(z), (_npad(dx), _npad(dy)), String(name))
end

Base.size(m::GravModel3D) = size(m.A)

# padding cells per side: those wider than the narrowest (core) cell
function _npad(d)
    core  = minimum(d)
    ncore = count(v -> v <= core * (1 + 1e-6), d)
    (length(d) - ncore) ÷ 2
end

"""
    sensor_depth(model, data) -> Vector

Depth (m) of each sensor below the model's depth zero, `top_elevation - elevation`; negative above the top of the mesh.
"""
sensor_depth(m::GravModel3D, d) = m.top_elevation .- d.elevation

"""
    air_mask(model) -> BitArray

Air cells (NaN) of a model.
"""
air_mask(m::GravModel3D) = isnan.(m.A)

"""
    cell_volumes(model) -> Array

Cell volumes (m³) in model order.
"""
cell_volumes(m::GravModel3D) = [a * b * c for a in m.dx, b in m.dy, c in m.dz]

"""
    surface_depth(model) -> Matrix

Depth (m) of the top of the first earth cell in each column, `Inf` for an all-air column.
"""
function surface_depth(m::GravModel3D)
    nx, ny, nz = size(m.A)
    s = fill(Inf, nx, ny)
    for i in 1:nx, j in 1:ny
        k = findfirst(!isnan, view(m.A, i, j, :))
        k === nothing || (s[i, j] = m.z[k])
    end
    s
end

#---------- reading ----------

"""
    ReadModel3D(path) -> GravModel3D

Read a density model written by `WriteModel3D` or `MakeMesh3D`.
"""
function ReadModel3D(path::AbstractString)
    isfile(path) || error("model file not found: $path")
    tokens = String[]
    header = ""
    top    = 0.0
    for raw in eachline(path)
        stripped = strip(raw)
        if startswith(stripped, '#')
            body = strip(stripped[2:end])
            if occursin(':', body)
                key, value = strip.(split(body, ':'; limit = 2))
                _header_key(key) == "top elevation (m)" && (top = parse(Float64, value))
            end
            continue
        end
        line = strip(first(split(raw, '#')))
        isempty(line) && continue
        if isempty(header)
            header = line
            continue
        end
        append!(tokens, split(line))
    end
    h = split(header)
    length(h) >= 3 && all(t -> tryparse(Int, t) !== nothing, h[1:3]) ||
        error("$path: expected a header 'nx ny nz nzAir LINEAR'")
    nx, ny, nz = parse.(Int, h[1:3])
    nzair = length(h) >= 4 ? something(tryparse(Int, h[4]), 0) : 0
    kind  = uppercase(length(h) >= 5 ? h[5] : (length(h) >= 4 && tryparse(Int, h[4]) === nothing ? h[4] : "LINEAR"))
    kind == "LINEAR" || error("$path: density models are LINEAR, got '$kind'")
    k = 0
    function next(n)
        k + n <= length(tokens) || error("$path: model file ends early")
        v = parse.(Float64, @view tokens[k+1:k+n])
        k += n
        v
    end
    dx, dy, dz = next(nx), next(ny), next(nz)
    A = Array{Float64, 3}(undef, nx, ny, nz)
    for iz in 1:nz
        A[:, :, iz] = reverse(reshape(next(nx * ny), nx, ny); dims = 1)
    end
    origin   = k + 3 <= length(tokens) ? next(3) : [0.0, 0.0, 0.0]
    rotation = k + 1 <= length(tokens) ? next(1)[1] : 0.0
    rotation == 0 || @warn "$path: rotation $rotation is ignored, the mesh is taken as north-aligned"
    A[A .> GRAV_AIR_THRESHOLD] .= NaN
    nzair == 0 || @warn "$path: $nzair air layers in the header are kept as model layers"
    GravModel3D(dx, dy, dz, A, origin; top_elevation = top, name = path)
end

#---------- writing ----------

"""
    WriteModel3D(path, model; title="") -> path

Write a model in the block layout, LINEAR density contrast (kg/m³), air as 1e17, with the top elevation as a header line.
"""
function WriteModel3D(path::AbstractString, m::GravModel3D; title::AbstractString = "")
    nx, ny, nz = size(m.A)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# Written by GravityGeophysics.jl WriteModel3D, density contrast (kg/m^3)", isempty(title) ? "" : ", $title")
        @printf(io, "# Top elevation (m) : %.6f\n", m.top_elevation)
        println(io, "$nx $ny $nz 0 LINEAR")
        for d in (m.dx, m.dy, m.dz)
            foreach(v -> print(io, v, " "), d)
            println(io)
        end
        for k in 1:nz
            println(io)
            for j in 1:ny
                for i in nx:-1:1
                    v = m.A[i, j, k]
                    @printf(io, "%15.6E", isnan(v) ? GRAV_AIR_VALUE : v)
                end
                println(io)
            end
        end
        println(io, m.origin[1], " ", m.origin[2], " ", m.origin[3])
        println(io, 0.0)
    end
    String(path)
end

"""
    model_like(model, values; name=model.name) -> GravModel3D

Copy of a model's mesh with new values (air kept as NaN).
"""
function model_like(m::GravModel3D, values::AbstractArray{<:Real, 3}; name::AbstractString = m.name)
    size(values) == size(m.A) || throw(DimensionMismatch("values are $(size(values)), the model $(size(m.A))"))
    A = Array{Float64, 3}(values)
    A[isnan.(m.A)] .= NaN
    GravModel3D(m.dx, m.dy, m.dz, A, m.origin; top_elevation = m.top_elevation, name)
end

"""
    same_mesh(a, b) -> Bool

True when two models share cell sizes, origin and top elevation.
"""
same_mesh(a::GravModel3D, b::GravModel3D) =
    size(a.A) == size(b.A) && a.dx ≈ b.dx && a.dy ≈ b.dy && a.dz ≈ b.dz && a.origin ≈ b.origin &&
    a.top_elevation ≈ b.top_elevation
