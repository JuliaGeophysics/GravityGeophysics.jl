# 3D gravity forward solver
# Author: @pankajkmishra
# Closed-form vertical attraction of uniform rectangular prisms (Plouff 1976; Blakely 1995; Nagy et al. 2000),
# evaluated once per mesh node and differenced over the eight corners of every cell, threaded over the data
# The sensitivity G is stored transposed, (n_cells × n_data) Float32, so each datum's row is contiguous
# Notation follows Blakely (1995): γ is the gravitational constant, ρ the density contrast, g the field
# A prism of density contrast ρ gives gz = -γ ρ Σ (-1)^(i+j+k) [ X log(Y + r) + Y log(X + r) - Z atan(XY/(Zr)) ],
# with X, Y, Z = corner - observation (x north, y east, z depth), r = |(X, Y, Z)| and the sum over the eight
# corners with the upper ones positive; gz is positive (downward) above a dense body. An atan term with Z = 0
# vanishes with its factor Z, and log(a + r) for a < 0 is evaluated as log((r² - a²)/(r - a)) to avoid cancellation

using LinearAlgebra
using JLD2
using Printf

"""
    GAMMA

γ, the gravitational constant (Blakely 1995), m³ kg⁻¹ s⁻².
"""
const GAMMA = 6.67430e-11

"""
    SI_TO_MGAL

Conversion from m/s² to mGal, 1e5.
"""
const SI_TO_MGAL = 1e5

#---------- kernel ----------

# log(a + r) without cancellation when a is large and negative; b2 = r² - a²
@inline _log_plus(a, b2, r) = a >= 0 ? log(a + r) : log(b2 / (r - a))
@inline _atan_ratio(num, den) = den == 0 ? 0.0 : atan(num / den)

"""
    GravKernel3D

Geometry shared by every datum of a survey on one mesh.
"""
struct GravKernel3D
    xn::Vector{Float64}
    yn::Vector{Float64}
    zn::Vector{Float64}
    obs::Matrix{Float64}
    scale::Float64
end

"""
    GravKernel3D(model, data)

Kernel of a survey on a model's mesh; sensor depths are `model.top_elevation - elevation`.
"""
function GravKernel3D(m::GravModel3D, d::GravData3D)
    nudge = 1e-6 .* (0.7071067811865476, 0.5773502691896258, 0.4472135954999579)
    obs   = permutedims(hcat(d.x .+ nudge[1], d.y .+ nudge[2], sensor_depth(m, d) .+ nudge[3]))
    GravKernel3D(m.x, m.y, m.z, obs, -GAMMA * SI_TO_MGAL)
end

ndata(k::GravKernel3D) = size(k.obs, 2)
ncells(k::GravKernel3D) = (length(k.xn) - 1) * (length(k.yn) - 1) * (length(k.zn) - 1)

# corner function at every node for datum i, into f (nx+1, ny+1, nz+1)
function _node_values!(f::Array{Float64, 3}, k::GravKernel3D, i::Int)
    ox, oy, oz = k.obs[1, i], k.obs[2, i], k.obs[3, i]
    @inbounds for (kk, zv) in enumerate(k.zn)
        Z  = zv - oz
        Z2 = Z * Z
        for (jj, yv) in enumerate(k.yn)
            Y  = yv - oy
            Y2 = Y * Y
            for (ii, xv) in enumerate(k.xn)
                X = xv - ox
                X2 = X * X
                r = sqrt(X2 + Y2 + Z2)
                f[ii, jj, kk] = X * _log_plus(Y, X2 + Z2, r) + Y * _log_plus(X, Y2 + Z2, r) - Z * _atan_ratio(X * Y, Z * r)
            end
        end
    end
    f
end

# corner sum of the node values for every cell, scaled to mGal per kg/m³, in model order
function _cell_row!(row::AbstractVector, f::Array{Float64, 3}, scale::Float64)
    nx, ny, nz = size(f) .- 1
    n = 0
    @inbounds for c in 1:nz, b in 1:ny, a in 1:nx
        n += 1
        s = f[a+1, b+1, c+1] - f[a, b+1, c+1] - f[a+1, b, c+1] + f[a, b, c+1] -
            f[a+1, b+1, c] + f[a, b+1, c] + f[a+1, b, c] - f[a, b, c]
        row[n] = scale * s
    end
    row
end

"""
    sensitivity_row(kernel, i) -> Vector{Float64}

Sensitivity of datum `i` to every cell, mGal per kg/m³, in model order.
"""
function sensitivity_row(k::GravKernel3D, i::Integer)
    f = Array{Float64, 3}(undef, length(k.xn), length(k.yn), length(k.zn))
    _cell_row!(Vector{Float64}(undef, ncells(k)), _node_values!(f, k, Int(i)), k.scale)
end

# run body(i, row) over the data, one node and row buffer per task
function _foreach_datum(body, k::GravKernel3D; progress::Bool = false, label::AbstractString = "")
    n      = ndata(k)
    nt     = max(1, min(Threads.nthreads(), n))
    done   = Threads.Atomic{Int}(0)
    step   = max(1, n ÷ 10)
    chunks = [i:nt:n for i in 1:nt]
    @sync for chunk in chunks
        Threads.@spawn begin
            f   = Array{Float64, 3}(undef, length(k.xn), length(k.yn), length(k.zn))
            row = Vector{Float64}(undef, ncells(k))
            for i in chunk
                _cell_row!(row, _node_values!(f, k, i), k.scale)
                body(i, row)
                c = Threads.atomic_add!(done, 1) + 1
                progress && c % step == 0 && (@printf("  %s %3d%%\n", label, round(Int, 100c / n)); flush(stdout))
            end
        end
    end
    nothing
end

#---------- forward ----------

"""
    forward3d(model, data, fwd::FwdCtrl3D=FwdCtrl3D()) -> Vector{Float64}

Predicted gz (mGal) of a model, air and NaN cells taken as zero, without storing G.
"""
function forward3d(m::GravModel3D, d::GravData3D, fwd::FwdCtrl3D = FwdCtrl3D())
    check_receivers(m, d)
    k = GravKernel3D(m, d)
    ρ = vec(replace(m.A, NaN => 0.0))
    pred = zeros(ndata(k))
    _foreach_datum((i, row) -> (pred[i] = dot(row, ρ)), k)
    pred
end

"""
    Sensitivity3D(model, data, fwd=FwdCtrl3D(); cells=nothing, storage=:memory, path="") -> Matrix{Float32}

Transposed sensitivity `Gᵀ`, `(n_cells × n_data)`, mGal per kg/m³, for the given linear cell indices (default every cell).
"""
function Sensitivity3D(m::GravModel3D, d::GravData3D, fwd::FwdCtrl3D = FwdCtrl3D();
                       cells::Union{Nothing, AbstractVector{<:Integer}} = nothing, verbose::Bool = false)
    check_receivers(m, d)
    k   = GravKernel3D(m, d)
    idx = cells === nothing ? collect(1:ncells(k)) : collect(Int, cells)
    nc  = length(idx)
    nd  = ndata(k)
    Gt  = Matrix{Float32}(undef, nc, nd)
    verbose && @printf("Sensitivity: %d data × %d cells, %.2f GB\n", nd, nc, 4.0 * nc * nd / 2^30)
    _foreach_datum(k; progress = verbose, label = "sensitivity") do i, row
        @inbounds for (j, c) in enumerate(idx)
            Gt[j, i] = row[c]
        end
    end
    Gt
end

"""
    SaveSensitivity3D(path, Gt; cells=nothing) -> path

Write the transposed sensitivity to a `.jld2` data file, with the shape and the linear cell indices it was built for.
"""
function SaveSensitivity3D(path::AbstractString, Gt::AbstractMatrix; cells = nothing)
    mkpath(dirname(abspath(path)))
    JLD2.jldsave(String(path); sensitivity = Matrix{Float32}(Gt), units = "mGal per kg/m3",
                 cells = cells === nothing ? collect(1:size(Gt, 1)) : collect(Int, cells))
    String(path)
end

"""
    LoadSensitivity3D(path) -> (sensitivity, cells, units)

Read a sensitivity file written by `SaveSensitivity3D`.
"""
function LoadSensitivity3D(path::AbstractString)
    isfile(path) || error("sensitivity file not found: $path")
    JLD2.jldopen(String(path), "r") do f
        (sensitivity = f["sensitivity"], cells = f["cells"], units = f["units"])
    end
end

#---------- receivers ----------

"""
    check_receivers(model, data)

Error when an observation point lies inside an earth cell; the prism formulas hold only outside the cells.
"""
function check_receivers(m::GravModel3D, d::GravData3D)
    nx, ny, nz = size(m.A)
    depth = sensor_depth(m, d)
    bad   = Int[]
    for i in eachindex(d.x)
        ix = searchsortedlast(m.x, d.x[i]); iy = searchsortedlast(m.y, d.y[i]); iz = searchsortedlast(m.z, depth[i] - 1e-6)
        (1 <= ix <= nx && 1 <= iy <= ny && 1 <= iz <= nz) || continue
        isnan(m.A[ix, iy, iz]) || push!(bad, i)
    end
    isempty(bad) && return nothing
    i = bad[1]
    error(@sprintf("%d observation(s) lie inside earth cells, first datum %d at x = %.1f, y = %.1f, elevation %.1f m (depth %.1f m below the mesh top); raise the sensors or rebuild the mesh with the topography",
                   length(bad), i, d.x[i], d.y[i], d.elevation[i], depth[i]))
end

#---------- file workflow ----------

"""
    ForwardSolve3D(model_path, data_path, fwd_path; output_path=<data dir>/data.pred) -> path

Predicted data of a model on the survey of a data file: the data file gives the locations and the errors, which are copied to the output.
"""
function ForwardSolve3D(model_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString;
                        output_path::AbstractString = joinpath(dirname(abspath(data_path)), "data.pred"))
    m, d, fwd = ReadModel3D(model_path), ReadData3D(data_path), ReadFwdCtrl3D(fwd_path)
    @printf("Forward: %d data on %d × %d × %d cells, %d threads\n", ndata(d), size(m.A)..., Threads.nthreads())
    elapsed = @elapsed pred = if fwd.write_frechet
        Gt = Sensitivity3D(m, d, fwd; verbose = true)
        SaveSensitivity3D(joinpath(dirname(abspath(output_path)), "sensitivity.jld2"), Gt)
        Float64.(Gt' * Float32.(vec(replace(m.A, NaN => 0.0))))
    else
        forward3d(m, d, fwd)
    end
    @printf("Forward: %.1f s\n", elapsed)
    WriteData3D(output_path, data_like(d, pred; title = "Predicted $(d.title)"))
end
