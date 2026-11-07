# 3D gravity regularized least-squares inversion by LSQR
# Author: @pankajkmishra
# Each step solves the regularized Gauss-Newton system of the Inv3D objective at the current β as the
# least-squares problem
#   min_δ ‖ A δ − b ‖²,   A = [ W_d J ;  √(β α_s w) ;  √(β α_i) W_i ],
#                          b = −[ W_d r ;  √(β α_s w) (x − x_ref) ;  √(β α_i) W_i x ]
# (W_d = 1/σ, r = F(ρ) − d), whose normal equations AᵀA δ = Aᵀb are (JᵀW_d²J + β H_m) δ = −∇Φ. LSQR (Paige &
# Saunders 1982) works on A and Aᵀ only, never forming AᵀA, so it is better conditioned than CG on the normal
# equations. Columns are scaled by the Jacobi factor 1/√diag(AᵀA) (right preconditioning), and only the cells that
# are off their bounds (or move inward) are free, so the driver's projected line search keeps the bounds.
# gz is linear in ρ, which makes each L2 step the exact least-squares solution at that β;
# the driver lowers β until the rms reaches the target (discrepancy principle). J is applied, never inverted, so
# the dense and matrix-free solvers take the same step.

using LinearAlgebra

#---------- configuration ----------

"""
    LSQR3DConfig(; max_iter=200, tolerance=1e-4)

Regularized least-squares algorithm for `Invert3D`, each step solved by LSQR.
- `max_iter`, `tolerance`: LSQR iterations of a solve, and its stop on ‖Aᵀr‖ / (‖A‖ ‖r‖)
"""
Base.@kwdef struct LSQR3DConfig <: AbstractInversion3D
    max_iter::Int = 200
    tolerance::Float64 = 1e-4
end

const LSQR3DResult = Inv3DResult{LSQR3DConfig}

inv3d_tag(::LSQR3DConfig) = "lsqr"

function inv3d_validate(alg::LSQR3DConfig)
    alg.max_iter > 0 || throw(ArgumentError("LSQR max_iter must be positive"))
    isfinite(alg.tolerance) && alg.tolerance > 0 || throw(ArgumentError("LSQR tolerance must be finite and positive"))
    nothing
end

#---------- interface ----------

mutable struct LSQR3DWork
    iterations::Int                             # LSQR iterations of the last direction
end

inv3d_init(alg::LSQR3DConfig, P, x) = LSQR3DWork(0)

inv3d_info(::LSQR3DConfig, work::LSQR3DWork) = (; lsqr = work.iterations)

inv3d_prepare!(alg::LSQR3DConfig, work::LSQR3DWork, P, x) = gradient(P, x)

# LSQR on min ‖A s − b‖ for s = D δ (column scaling D), A given by its products; returns δ and the iterations
function _lsqr(Av, Atu, b, D, max_iter, tol)
    n = length(D)
    u = copy(b); β = norm(u)
    s = zeros(n)
    β == 0 && return s, 0
    u ./= β
    v = D .* Atu(u); α = norm(v)
    α == 0 && return s, 0
    v ./= α
    w = copy(v)
    φbar, ρbar = β, α
    normA2 = α^2
    k = 0
    while k < max_iter
        k += 1
        u = Av(D .* v) .- α .* u; β = norm(u)
        β > 0 && (u ./= β)
        normA2 += β^2
        v = D .* Atu(u) .- β .* v; α = norm(v)
        α > 0 && (v ./= α)
        normA2 += α^2
        ρ = hypot(ρbar, β)
        c, sn = ρbar / ρ, β / ρ
        θ = sn * α
        ρbar = -c * α
        φ = c * φbar
        φbar = sn * φbar
        s .+= (φ / ρ) .* w
        w = v .- (θ / ρ) .* w
        # ‖Aᵀr‖ = φbar α |c|, ‖r‖ = φbar
        (α == 0 || φbar * α * abs(c) <= tol * sqrt(normA2) * φbar) && break
    end
    D .* s, k
end

# the stacked system at x over the free cells: products with A and Aᵀ, and b
function _lsqr3d_system(P, x, free)
    ρ = _full(P, x)
    sb = sqrt(P.beta)
    ws = sb .* sqrt.(P.alpha_s .* P.w)
    ri = [sb * sqrt(P.alpha[i]) for i in 1:3]
    nd = length(P.d)
    off = cumsum([nd, length(x), (size(Wi, 1) for Wi in P.W)...])  # ends of the blocks: data, smallness, x, y, z
    function Av(v)
        f = free .* v
        vcat(jvec(P.op, ρ, _embed(P, f)) ./ P.σ, ws .* f, (ri[i] .* (P.W[i] * f) for i in 1:3)...)
    end
    function Atu(u)
        t = jtvec(P.op, ρ, u[1:nd] ./ P.σ)[P.var] .+ ws .* u[nd+1:off[2]]
        for i in 1:3
            t .+= P.W[i]' * (ri[i] .* u[off[i+1]+1:off[i+2]])
        end
        free .* t
    end
    b = -vcat((forward(P.op, ρ) .- P.d) ./ P.σ, ws .* (x .- P.xref), (ri[i] .* (P.W[i] * x) for i in 1:3)...)
    Av, Atu, b
end

# one solve at x over the free cells
function _lsqr3d_solve(alg::LSQR3DConfig, P, x, g)
    free = .!(((x .<= P.lower) .& (g .> 0)) .| ((x .>= P.upper) .& (g .< 0)))
    Av, Atu, b = _lsqr3d_system(P, x, free)
    # Jacobi column scaling: diag(AᵀA) = diag JᵀJ/σ² (at the start) + β(α_s w + Σ α_i diag W_iᵀW_i)
    d = P.jtj .+ P.beta .* (P.alpha_s .* P.w .+ sum(P.alpha[i] .* P.WtWdiag[i] for i in 1:3))
    D = free ./ sqrt.(max.(d, eps()))
    _lsqr(Av, Atu, b, D, alg.max_iter, alg.tolerance)
end

function inv3d_direction(alg::LSQR3DConfig, work::LSQR3DWork, P, x, g)
    step, work.iterations = _lsqr3d_solve(alg, P, x, g)
    step
end

#---------- entry points ----------

"""
    LSQR3D(m, d; config=LSQR3DConfig(), options=Inv3DOptions(), fwd=FwdCtrl3D(), mask=trues, reference=nothing)

LSQR shorthand for `Invert3D(...; algorithm=config)`.
"""
LSQR3D(m::Model, d::Data; config::LSQR3DConfig = LSQR3DConfig(), kwargs...) =
    Invert3D(m, d; algorithm = config, kwargs...)
