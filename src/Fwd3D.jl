# 3D gravity forward modelling with prisms
# Author: @pankajkmishra
# Closed-form vertical attraction of uniform right rectangular prisms (Plouff 1976; Blakely 1995; Nagy et al. 2000).
# Notation follows Blakely (1995): γ is the gravitational constant, ρ the density contrast. A prism gives
#   gz = −γ ρ Σ (−1)^(i+j+k) [ X log(Y + r) + Y log(X + r) − Z atan(XY / (Z r)) ],
# X, Y, Z = vertex − station (x north, y east, z DOWN), r = |(X, Y, Z)|, the sum over the eight vertices with the
# upper ones positive; gz is positive (down) above a dense body. The vertex function depends on the vertex only,
# so it is evaluated once per mesh node and differenced over the eight corners of every cell. An atan term with
# Z = 0 vanishes with its factor Z, and log(a + r) for a < 0 is evaluated as log((r² − a²)/(r − a)) against
# cancellation. On a mesh plane a term is 0 · log 0 and is taken as its limit 0 (x log x → 0), so stations may lie
# on cell faces and edges, as ground stations on the earth cells do.
# One row routine serves both solvers, so the dense G and the matrix-free products are identical.

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

#---------- kernel (x north, y east, z down; vertex − station) ----------

# log(a + r) without cancellation when a is large and negative; b2 = r² − a²
@inline _log_plus(a, b2, r) = a >= 0 ? log(a + r) : log(b2 / (r - a))
# c log(a + r), 0 when the factor c is 0 (the limit on a mesh plane)
@inline _xlog(c, a, b2, r) = c == 0 ? 0.0 : c * _log_plus(a, b2, r)
@inline _atan_ratio(num, den) = den == 0 ? 0.0 : atan(num / den)

# S[node] = vertex function of station (x0, y0, z0)
function _node_values!(S::Array{Float64,3}, m::Model, x0, y0, z0)
    @inbounds for k in 1:m.nz+1
        Z = m.z[k] - z0
        Z2 = Z * Z
        for j in 1:m.ny+1
            Y = m.y[j] - y0
            Y2 = Y * Y
            for i in 1:m.nx+1
                X = m.x[i] - x0
                X2 = X * X
                r = sqrt(X2 + Y2 + Z2)
                S[i, j, k] = _xlog(X, Y, X2 + Z2, r) + _xlog(Y, X, Y2 + Z2, r) - Z * _atan_ratio(X * Y, Z * r)
            end
        end
    end
    S
end

# row of the active cells: −γ × 1e5 Σ over the 8 vertices, upper positive; mGal per kg/m³
function _row!(row::AbstractVector, S::Array{Float64,3}, cells::Vector{CartesianIndex{3}})
    c = -GAMMA * SI_TO_MGAL
    @inbounds for (n, I) in enumerate(cells)
        i, j, k = I[1], I[2], I[3]
        row[n] = c * (S[i+1, j+1, k+1] - S[i, j+1, k+1] - S[i+1, j, k+1] + S[i, j, k+1] -
                      S[i+1, j+1, k] + S[i, j+1, k] + S[i+1, j, k] - S[i, j, k])
    end
    row
end

#---------- operator ----------

"""
    GravForward(m, d; solver=:dense)

Forward operator of the data `d` on the active (non-air) cells of `m`. `solver = :dense` builds and
stores the sensitivity G (mGal per kg/m³); `:matrixfree` recomputes the rows for every product and
needs no storage. Both give the same numbers. Checks the data frame against the mesh.
"""
struct GravForward
    m::Model
    d::Data
    cells::Vector{CartesianIndex{3}}
    zs::Vector{Float64}                   # station z on the mesh (d.z + m.datum)
    solver::Symbol
    G::Matrix{Float64}                    # 0 × 0 for :matrixfree
end

function GravForward(m::Model, d::Data; solver::Symbol = :dense)
    solver in (:dense, :matrixfree) || throw(ArgumentError("solver must be :dense or :matrixfree"))
    check_frame(d, m)
    cells = findall(active_cells(m))
    zs = _mesh_z(d, m)
    solver == :matrixfree && return GravForward(m, d, cells, zs, solver, zeros(0, 0))
    bytes = 8.0 * ndata(d) * length(cells)
    bytes > 0.8 * Sys.total_memory() &&
        throw(ArgumentError(@sprintf("dense G needs %.1f GB; use solver = :matrixfree", bytes / 2^30)))
    op = GravForward(m, d, cells, zs, solver, Matrix{Float64}(undef, ndata(d), length(cells)))
    _foreach_row(op) do n, row
        op.G[n, :] .= row
    end
    op
end

ncells(op::GravForward) = length(op.cells)

# run f(n, row) for every station n, threaded over chunks of stations
function _foreach_row(f, op::GravForward)
    m, d = op.m, op.d
    N = ndata(d)
    chunks = collect(Iterators.partition(1:N, cld(N, max(1, Threads.nthreads()))))
    @sync for chunk in chunks
        Threads.@spawn begin
            S = Array{Float64,3}(undef, m.nx + 1, m.ny + 1, m.nz + 1)
            row = Vector{Float64}(undef, length(op.cells))
            for n in chunk
                _node_values!(S, m, d.x[n], d.y[n], op.zs[n])
                f(n, _row!(row, S, op.cells))
            end
        end
    end
    nothing
end

# Σ_n f(a, n, row) into one accumulator per chunk of stations, then summed
function _accumulate_rows(f, op::GravForward)
    acc = [zeros(ncells(op)) for _ in 1:max(1, Threads.nthreads())]
    m, d = op.m, op.d
    N = ndata(d)
    chunks = collect(Iterators.partition(1:N, cld(N, length(acc))))
    @sync for (t, chunk) in enumerate(chunks)
        Threads.@spawn begin
            S = Array{Float64,3}(undef, m.nx + 1, m.ny + 1, m.nz + 1)
            row = Vector{Float64}(undef, ncells(op))
            a = acc[t]
            for n in chunk
                _node_values!(S, m, d.x[n], d.y[n], op.zs[n])
                f(a, n, _row!(row, S, op.cells))
            end
        end
    end
    sum(acc)
end

# G v
function _gv(op::GravForward, v::AbstractVector)
    op.solver == :dense && return op.G * v
    out = zeros(ndata(op.d))
    _foreach_row(op) do n, row
        out[n] = dot(row, v)
    end
    out
end

# Gᵀ w
function _gtw(op::GravForward, w::AbstractVector)
    op.solver == :dense && return transpose(op.G) * w
    _accumulate_rows((a, n, row) -> (a .+= w[n] .* row), op)
end

_check_len(op, v) = length(v) == ncells(op) ||
    throw(ArgumentError("model vector has $(length(v)) values, the operator $(ncells(op)) active cells"))

"""
    forward(op, ρ) -> predicted data

`ρ` holds the density contrast (kg/m³) of the active cells, in `findall(active_cells(m))` order;
the result is gz in mGal.
"""
forward(op::GravForward, ρ::AbstractVector) = (_check_len(op, ρ); _gv(op, ρ))

"""
    jvec(op, ρ, v) / jtvec(op, ρ, w)

Jacobian products at ρ. gz is linear in ρ: J = G whatever ρ.
"""
jvec(op::GravForward, ρ::AbstractVector, v::AbstractVector) = (_check_len(op, v); _gv(op, v))

function jtvec(op::GravForward, ρ::AbstractVector, w::AbstractVector)
    length(w) == ndata(op.d) || throw(ArgumentError("data vector length differs from the data"))
    _gtw(op, w)
end

"""
    model_vector(op, m) / model_array(op, ρ)

Active-cell vector of a model on the operator's mesh, and the full array (NaN in air) of a vector.
"""
model_vector(op::GravForward, m::Model) = (same_mesh(op.m, m) || throw(ArgumentError("model mesh differs")); m.A[op.cells])
function model_array(op::GravForward, ρ::AbstractVector)
    _check_len(op, ρ)
    A = fill(NaN, size(op.m.A)); A[op.cells] .= ρ; A
end

#---------- file workflow ----------

"""
    ForwardSolve3D(model_path, data_path, fwd_path; output_path, add_noise=false, rng_seed=20261017) -> path

Forward run from files: a density-contrast model, a native data file giving the stations and errors,
and `FwdCtrl` giving the solver. Writes the predicted data, by default `data.pred` next to the data
file, with the errors of the data file; `add_noise` adds Gaussian noise of those errors.
"""
function ForwardSolve3D(model_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString;
                        output_path::AbstractString = joinpath(dirname(abspath(data_path)), "data.pred"),
                        add_noise::Bool = false, rng_seed::Integer = 20261017)
    fwd = ReadFwdCtrl3D(fwd_path)
    m, d = load_model(model_path), load_data_grav(data_path)
    op = GravForward(m, d; solver = fwd.solver)
    p = forward(op, model_vector(op, m))
    add_noise && (p .+= d.error .* randn(Xoshiro(rng_seed), length(p)))
    write_data_grav(output_path, d; value = p, comment = "predicted by $(basename(model_path))" * (add_noise ? ", with noise" : ""))
end
