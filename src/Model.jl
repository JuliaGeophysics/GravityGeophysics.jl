# Model file reader/writer and the shared gravity mesh.
# Author: @pankajkmishra
# x north, y east, z DOWN. On file: header, cell widths, one block per layer with x stored
# backwards within each y row, then the origin (x, y, z of the first corner) and the rotation.
# Vertical datum: z = 0 lies at `datum` m above sea level, written as the comment line
# "# Datum elevation (m) : value". MakeMesh3D puts z = 0 at the highest ground (the top of the highest
# earth cell); a file without the line (MT, magnetic, older models) has its z = 0 at sea level.
# Data files keep z = −elevation (sea level); a station's z on the mesh is d.z + datum.
# Values are LINEAR density contrast in kg/m³; air and inactive cells are NaN in memory and AIR_VALUE
# on file, read back as NaN above AIR_THRESHOLD, the tag the resistivity files use for air.

#---------- model ----------

const AIR_VALUE = 1e17        # written for NaN cells
const AIR_THRESHOLD = 1e15    # read back as NaN above this

# Colour map of every plot: Spectral reversed, blue low to red high
const COLORMAP = Reverse(:Spectral)

mutable struct Model
    A::Array{Float64,3}       # density contrast (kg/m³), NaN = air
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    nx::Int
    ny::Int
    nz::Int
    x::Vector{Float64}        # edges
    y::Vector{Float64}
    z::Vector{Float64}
    origin::Vector{Float64}
    cx::Vector{Float64}       # centres
    cy::Vector{Float64}
    cz::Vector{Float64}
    X::Matrix{Float64}
    Y::Matrix{Float64}
    Xc::Matrix{Float64}
    Yc::Matrix{Float64}
    Z::Matrix{Float64}        # top of the earth in each column
    npad::NTuple{2,Int}
    name::String
    niter::String
    rotation::Float64         # degrees, mesh x clockwise from north
    datum::Float64            # elevation (m a.s.l.) of z = 0
end


meshgrid(ax::AbstractVector, ay::AbstractVector) = (repeat(reshape(ax, 1, :), length(ay), 1), repeat(reshape(ay, :, 1), 1, length(ax)))

"""
    Model(dx, dy, dz, A=zeros; origin=[-Σdx/2, -Σdy/2, 0], rotation=0, datum=0, name="")

Density-contrast model on a tensor mesh, with its derived edge, centre and surface arrays.
`A[nx, ny, nz]` is density contrast in kg/m³, NaN in air. `datum` is the elevation (m a.s.l.) of z = 0.
"""
function Model(dx::AbstractVector, dy::AbstractVector, dz::AbstractVector,
               A::AbstractArray{<:Real,3} = zeros(length(dx), length(dy), length(dz));
               origin::AbstractVector = [-sum(dx)/2, -sum(dy)/2, 0.0], rotation::Real = 0.0, datum::Real = 0.0,
               name::AbstractString = "", niter::AbstractString = "")
    nx, ny, nz = length(dx), length(dy), length(dz)
    all(>(0), dx) && all(>(0), dy) && all(>(0), dz) || throw(ArgumentError("cell widths must be positive"))
    size(A) == (nx, ny, nz) || throw(ArgumentError("model size $(size(A)) differs from mesh ($nx, $ny, $nz)"))
    length(origin) == 3 || throw(ArgumentError("origin must have 3 values"))
    x = vcat(0.0, cumsum(dx)) .+ origin[1]
    y = vcat(0.0, cumsum(dy)) .+ origin[2]
    z = vcat(0.0, cumsum(dz)) .+ origin[3]
    c(e) = (e[1:end-1] .+ e[2:end]) ./ 2
    cx, cy, cz = c(x), c(y), c(z)
    X, Y = meshgrid(y, x)
    Xc, Yc = meshgrid(cy, cx)
    m = Model(Array{Float64}(A), collect(Float64, dx), collect(Float64, dy), collect(Float64, dz), nx, ny, nz,
              x, y, z, collect(Float64, origin), cx, cy, cz, X, Y, Xc, Yc, zeros(nx, ny), (0, 0),
              String(name), String(niter), Float64(rotation), Float64(datum))
    m.Z = _surface(m)
    m.npad = _npad(m)
    m
end

# z of the first earth cell under the top run of air cells; a column of air only keeps the
# top of its deepest cell
function _surface(m::Model)
    Z = zeros(m.nx, m.ny)
    @inbounds for i in 1:m.nx, j in 1:m.ny
        ind = something(findlast(isnan, view(m.A, i, j, :)), 0)
        Z[i, j] = m.z[min(ind, m.nz - 1) + 1]
    end
    Z
end

# padding cells per side: cells whose width differs from the central one
function _npad(m::Model)
    cx0 = m.dx[clamp(round(Int, m.nx/2), 1, m.nx)]
    cy0 = m.dy[clamp(round(Int, m.ny/2), 1, m.ny)]
    ((m.nx - count(==(cx0), m.dx)) ÷ 2, (m.ny - count(==(cy0), m.dy)) ÷ 2)
end

air_cells(m::Model) = isnan.(m.A)
active_cells(m::Model) = .!isnan.(m.A)

function same_mesh(a::Model, b::Model)
    a.dx == b.dx && a.dy == b.dy && a.dz == b.dz && a.origin == b.origin && a.rotation == b.rotation && a.datum == b.datum
end

#---------- model files ----------

# the datum comment line, 0 (sea level) when absent
function _datum(raw)
    for l in raw
        m = match(r"^\s*#\s*datum elevation \(m\)\s*:\s*(\S+)"i, l)
        m === nothing || return parse(Float64, m[1])
    end
    0.0
end

function _model_read(path::AbstractString)
    raw = readlines(path)
    datum = _datum(raw)
    lines = [strip(first(split(s, '#'))) for s in raw]
    filter!(!isempty, lines)
    isempty(lines) && throw(ArgumentError("empty model file: $path"))
    h = split(popfirst!(lines))
    length(h) == 5 || throw(ArgumentError("expected header 'nx ny nz nzAir LINEAR/LOGE/LOG10' in $path"))
    any(occursin("VALUES", uppercase(t)) for t in h) &&
        throw(ArgumentError("VALUES-format models are not supported: $path"))
    nx, ny, nz, nair = parse.(Int, h[1:4])
    all(>(0), (nx, ny, nz)) && 0 <= nair < nz || throw(ArgumentError("invalid model dimensions in $path"))
    scale = uppercase(h[5])
    scale in ("LINEAR", "LOGE", "LOG10") || throw(ArgumentError("unsupported model scale $scale"))
    tokens = split(join(lines, ' '))
    n = nx + ny + nz + nx*ny*nz
    length(tokens) in (n, n+3, n+4) ||
        throw(ArgumentError("expected $n values (+ origin, rotation) in $path, found $(length(tokens))"))
    v = parse.(Float64, tokens)
    all(isfinite, v) || throw(ArgumentError("nonfinite model values in $path"))
    dx, dy, dz = v[1:nx], v[nx+1:nx+ny], v[nx+ny+1:nx+ny+nz]
    origin = length(v) >= n+3 ? v[n+1:n+3] : [-sum(dx)/2, -sum(dy)/2, 0.0]
    rotation = length(v) == n+4 ? v[end] : 0.0
    A = reverse(reshape(v[nx+ny+nz+1:n], nx, ny, nz); dims = 1)
    linear = scale == "LINEAR" ? A : scale == "LOGE" ? exp.(A) : exp10.(A)
    linear[linear .> AIR_THRESHOLD] .= NaN
    (; dx, dy, dz, A = linear, origin, rotation, datum, nair, scale)
end

"""
    load_model(path) -> Model

Read a density model file (LINEAR, kg/m³). Cells above 1e15 are air (NaN).
Use [`load_mesh`](@ref) to start from a resistivity model.
"""
function load_model(path::AbstractString)
    r = _model_read(path)
    r.scale == "LINEAR" ||
        throw(ArgumentError("density models are LINEAR; $path is $(r.scale), use load_mesh for a resistivity model"))
    Model(r.dx, r.dy, r.dz, r.A; origin = r.origin, rotation = r.rotation, datum = r.datum, name = path)
end

"""
    load_mesh(path; value=0.0, water_log10=NaN) -> Model

Mesh, air cells and rotation of any model file (e.g. a resistivity model with topography),
filled with density contrast `value` below the surface. Air (> 1e15 Ω·m) stays NaN. Cells with
log10 ρ < `water_log10` (NaN: off) are set to NaN too: sea water carries no modelled contrast and must not
be inverted.
"""
function load_mesh(path::AbstractString; value::Real = 0.0, water_log10::Real = NaN)
    r = _model_read(path)
    A = map(r.A) do rho
        isnan(rho) || (!isnan(water_log10) && log10(rho) < water_log10) ? NaN : Float64(value)
    end
    Model(r.dx, r.dy, r.dz, A; origin = r.origin, rotation = r.rotation, datum = r.datum, name = path)
end

# top layers made of air only
_top_air_layers(A) = something(findfirst(k -> !all(isnan, view(A, :, :, k)), axes(A, 3)), size(A, 3)) - 1

"""
    write_model(path, m; nzAir=<air-only top layers>, comment="")

Write a density model file: header, cell widths, one block per layer with x running
backwards within each y row, then origin and rotation. Values are LINEAR kg/m³; NaN cells are
written as $(AIR_VALUE). The datum is written as the line "# Datum elevation (m) : value".
"""
function write_model(path::AbstractString, m::Model; nzAir::Integer = _top_air_layers(m.A),
                           comment::AbstractString = "")
    size(m.A) == (m.nx, m.ny, m.nz) || throw(ArgumentError("model shape differs from mesh"))
    all(v -> isnan(v) || isfinite(v), m.A) || throw(ArgumentError("density contrast must be finite or NaN (air)"))
    0 <= nzAir < m.nz || throw(ArgumentError("invalid air layer count $nzAir"))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# Written by GravityGeophysics.jl write_model (density contrast, kg/m³)",
                isempty(comment) ? "" : " -- " * comment)
        println(io, "# Datum elevation (m) : ", m.datum)
        @printf(io, "%d %d %d %d LINEAR\n", m.nx, m.ny, m.nz, nzAir)
        for w in (m.dx, m.dy, m.dz)
            for v in w; print(io, v, " "); end; println(io)
        end
        for k in 1:m.nz
            println(io)
            for j in 1:m.ny
                for i in m.nx:-1:1
                    v = m.A[i, j, k]
                    @printf(io, "%17.9E", isnan(v) ? AIR_VALUE : v)
                end
                println(io)
            end
        end
        println(io, m.origin[1], " ", m.origin[2], " ", m.origin[3])
        println(io, m.rotation)
    end
    String(path)
end
