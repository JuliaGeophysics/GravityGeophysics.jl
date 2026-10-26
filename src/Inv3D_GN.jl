# 3D gravity Gauss–Newton inversion
# Author: @pankajkmishra
# Projected, damped (Levenberg–Marquardt) Gauss–Newton direction for the Invert3D driver
# The step solves (JᵀJ + βRᵀR + λD) p = -g on the cells off the bounds by Jacobi-preconditioned CG, matrix free

using LinearAlgebra

#---------- configuration ----------

"""
    GaussNewton3DConfig(; damping=1e-3, max_damping_trials=6, cg_iterations=50, cg_tolerance=1e-3)

Gauss–Newton algorithm for `Invert3D`.
"""
Base.@kwdef struct GaussNewton3DConfig <: AbstractInversion3D
    damping::Float64 = 1e-3
    max_damping_trials::Int = 6
    cg_iterations::Int = 50
    cg_tolerance::Float64 = 1e-3
end

const GaussNewton3DResult = Inv3DResult{GaussNewton3DConfig}

inv3d_tag(::GaussNewton3DConfig) = "gn"

function inv3d_validate(alg::GaussNewton3DConfig)
    isfinite(alg.damping) && alg.damping >= 0 || throw(ArgumentError("damping must be finite and nonnegative"))
    alg.max_damping_trials > 0 || throw(ArgumentError("max_damping_trials must be positive"))
    alg.cg_iterations > 0 && alg.cg_tolerance > 0 || throw(ArgumentError("CG iterations and tolerance must be positive"))
    nothing
end

#---------- interface ----------

mutable struct GaussNewton3DWork
    damping::Float64
    trials::Int
    JtJ::Vector{Float64}      # diag(JᵀJ), fixed
    diagonal::Vector{Float64} # diag(JᵀJ + βRᵀR) of the current β
    cg::Int                   # CG iterations of the last direction
end

inv3d_init(alg::GaussNewton3DConfig, problem, state) =
    GaussNewton3DWork(alg.damping, 0, _inv3d_JtJ_diag(problem.Gt, problem.sigma), Float64[], 0)

inv3d_info(::GaussNewton3DConfig, work::GaussNewton3DWork) = (; damping = work.damping, cg = work.cg)

function inv3d_prepare!(alg::GaussNewton3DConfig, work::GaussNewton3DWork, problem, state)
    work.diagonal = work.JtJ .+ state.beta .* diag(problem.RtR)
    work.trials   = 0
    inv3d_gradient(problem, state)
end

function inv3d_direction(alg::GaussNewton3DConfig, work::GaussNewton3DWork, problem, state, gradient)
    free = .!inv3d_bound_cells(problem, state, gradient)
    d    = max.(work.diagonal, 1e-30)
    mul(v) = (inv3d_hessian_mul(problem, state, v) .+ work.damping .* d .* v) .* free
    p, work.cg = _pcg(mul, -gradient .* free, free ./ ((1 + work.damping) .* d); maxiter = alg.cg_iterations,
                      rtol = alg.cg_tolerance)
    p
end

function inv3d_reject!(alg::GaussNewton3DConfig, work::GaussNewton3DWork)
    work.damping = max(work.damping, 1e-4) * 10
    work.trials += 1
    work.trials < alg.max_damping_trials
end

function inv3d_accept!(alg::GaussNewton3DConfig, work::GaussNewton3DWork, problem, old, new)
    work.damping = work.damping / 3
    nothing
end

#---------- conjugate gradients ----------

# preconditioned CG for A x = b, A given as a function and Minv the inverse diagonal
function _pcg(A, b, Minv; maxiter::Int, rtol::Real)
    x     = zeros(length(b))
    r     = copy(b)
    bnorm = norm(b)
    bnorm == 0 && return x, 0
    z  = Minv .* r
    p  = copy(z)
    rz = dot(r, z)
    for it in 1:maxiter
        Ap  = A(p)
        pAp = dot(p, Ap)
        pAp <= 0 && return x, it
        alpha = rz / pAp
        x .+= alpha .* p
        r .-= alpha .* Ap
        norm(r) <= rtol * bnorm && return x, it
        z      = Minv .* r
        rz_new = dot(r, z)
        p .= z .+ (rz_new / rz) .* p
        rz = rz_new
    end
    x, maxiter
end

#---------- entry points ----------

"""
    GaussNewton3D(model, observed, fwd; config=GaussNewton3DConfig(), kwargs...)

Gauss–Newton shorthand for `Invert3D(...; algorithm=config)`.
"""
GaussNewton3D(model::GravModel3D, observed::GravData3D, fwd::FwdCtrl3D;
              config::GaussNewton3DConfig = GaussNewton3DConfig(), kwargs...) =
    Invert3D(model, observed, fwd; algorithm = config, kwargs...)
