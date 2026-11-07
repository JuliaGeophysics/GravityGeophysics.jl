# 3D gravity deterministic inversion framework
# Author: @pankajkmishra
# Algorithm-agnostic problem setup, objective, regularization, line search, and driver loop
# Each algorithm lives in its own Inv3D_<Name>.jl and plugs in through the interface below
#
# Unknowns x: density contrast of the variable cells (active and in the mask); fixed cells keep the
# starting values. The objective is
#   Φ = φ_d + β φ_m,   φ_d = ½ ‖(F(ρ) − d) / σ‖²,
#   φ_m = ½ α_s ‖√w ⊙ (x − x_ref)‖² + ½ Σ_i α_i ‖W_i x‖²
# (L2 smallness and smoothness). W_i are first differences between neighbouring variable cells along x, y, z, scaled by
# min(h)/h of the face and by √(mean w of the two cells). w are the cell weights: sensitivity
# (√ diag JᵀJ/σ², normalised), Li & Oldenburg depth weighting (depth below the stations)^(−p), or 1.
# β starts from the Rayleigh quotient ratio of the data and model Hessians times `beta` (the control
# file's lambda) and is divided by the cooling factor until the rms reaches the target; a step that
# overshoots it by more than the tolerance is redone with β raised, then bisected in log β, so the run
# stops at the target instead of fitting noise. gz is linear in the density, so J = G throughout.

using LinearAlgebra
using Printf
using SparseArrays
using Statistics

#---------- algorithm interface ----------
#
# an algorithm is a config struct `A <: AbstractInversion3D` plus these methods:
#
#   inv3d_init(alg, problem, x)                        -> work   (mutable scratch)
#   inv3d_prepare!(alg, work, problem, x)              -> gradient, once per step
#   inv3d_direction(alg, work, problem, x, grad)       -> search direction in density
#   inv3d_reject!(alg, work)                           -> true = retry with a new direction
#   inv3d_accept!(alg, work, problem, old, new)        -> nothing
#   inv3d_info(alg, work)                              -> NamedTuple added to history
#   inv3d_tag(alg)                                     -> short label for the log, e.g. "lsqr"
#   inv3d_validate(alg)                                -> throws on bad settings
#
# the driver owns stopping tests, step capping, the projected Armijo line search and the
# cooling of β, so a new algorithm only decides which direction to take

"""
    AbstractInversion3D

Supertype of 3D deterministic inversion algorithms. See `Inv3D.jl` for the interface.
"""
abstract type AbstractInversion3D end

inv3d_reject!(::AbstractInversion3D, work) = false
inv3d_accept!(::AbstractInversion3D, work, problem, old, new) = nothing
inv3d_info(::AbstractInversion3D, work) = (;)
inv3d_validate(::AbstractInversion3D) = nothing

#---------- options and result ----------

"""
    Inv3DOptions(; kwargs...)

Algorithm-independent inversion controls.
- `max_iter`: maximum accepted iterations; 0 evaluates the start model only
- `beta`: initial regularization weight as a multiple of the data/model Hessian ratio
- `cooling`, `cooling_interval`: β is divided by `cooling` every `cooling_interval` iterations
- `target_rms`, `target_tolerance`: stop at the target rms; a step below
  `(1 − target_tolerance) × target_rms` is redone with a larger β
- `smallness`, `smooth_x`, `smooth_y`, `smooth_z`: weights of the regularization terms
- `weighting`, `depth_exponent`: cell weights `:depth` (default, exponent 2), `:sensitivity` or `:none`
- `bounds`: density bounds (kg/m³), kept by projection; per-cell bounds (`Invert3D` `lower`, `upper`) replace them
- `max_step`: cap on the largest density change per iteration (kg/m³)
- `max_linesearch`: backtracking halvings per direction
- `verbose`: print one line per accepted iteration
"""
Base.@kwdef struct Inv3DOptions
    max_iter::Int = 30
    beta::Float64 = 1.0
    cooling::Float64 = 2.0
    cooling_interval::Int = 1
    target_rms::Float64 = 1.0
    target_tolerance::Float64 = 0.1
    smallness::Float64 = 0.01
    smooth_x::Float64 = 1.0
    smooth_y::Float64 = 1.0
    smooth_z::Float64 = 1.0
    weighting::Symbol = :depth
    depth_exponent::Float64 = 2.0
    bounds::Tuple{Float64, Float64} = (-1000.0, 1000.0)
    max_step::Float64 = Inf
    max_linesearch::Int = 20
    verbose::Bool = true
end

"""
    Inv3DResult

Result of `Invert3D`.
- `model`: recovered density-contrast model (kg/m³), air included
- `predicted`, `rms`: predicted data and their error-weighted rms
- `history`: accepted iterations, iteration zero first
- `active_cells`: inverted cells (indices into the model array)
- `converged`: false for `:max_iter` and `:line_search_failed`
- `reason`: termination reason
- `algorithm`: the algorithm config that produced it
- `problem`: the `Inv3DProblem`, for the objective functions and joint inversion
"""
struct Inv3DResult{A<:AbstractInversion3D}
    model::Model
    predicted::Vector{Float64}
    rms::Float64
    history::Vector{NamedTuple}
    active_cells::Vector{CartesianIndex{3}}
    converged::Bool
    reason::Symbol
    algorithm::A
    problem::Any
end

function _inv3d_validate(opt::Inv3DOptions)
    opt.max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    opt.weighting in (:none, :depth, :sensitivity) || throw(ArgumentError("weighting must be :none, :depth or :sensitivity"))
    all(x -> isfinite(x) && x >= 0, (opt.smallness, opt.smooth_x, opt.smooth_y, opt.smooth_z, opt.target_rms,
                                     opt.depth_exponent)) ||
        throw(ArgumentError("regularization and stopping controls must be finite and nonnegative"))
    isfinite(opt.beta) && opt.beta > 0 || throw(ArgumentError("beta must be finite and positive"))
    opt.cooling > 1 || throw(ArgumentError("the cooling factor must exceed 1"))
    opt.cooling_interval > 0 || throw(ArgumentError("cooling_interval must be positive"))
    0 <= opt.target_tolerance < 1 || throw(ArgumentError("target_tolerance must lie in [0, 1)"))
    all(isfinite, opt.bounds) && opt.bounds[1] < opt.bounds[2] || throw(ArgumentError("invalid density bounds"))
    opt.max_step > 0 || throw(ArgumentError("max_step must be positive"))
    opt.max_linesearch > 0 || throw(ArgumentError("max_linesearch must be positive"))
    nothing
end

#---------- problem ----------

"""
    Inv3DProblem

Fixed part of an inversion: forward operator, data, variable cells, regularization operators,
reference model, bounds and the current β. Built once by `Invert3D`.
"""
mutable struct Inv3DProblem
    op::GravForward
    d::Vector{Float64}
    σ::Vector{Float64}
    var::Vector{Int}               # operator columns of the variable cells
    ρ::Vector{Float64}             # all active cells; fixed ones stay at their start values
    xref::Vector{Float64}
    w::Vector{Float64}             # cell weights on the variable cells
    W::Vector{SparseMatrixCSC{Float64,Int}}
    WtWdiag::Vector{Vector{Float64}}
    alpha_s::Float64
    alpha::NTuple{3,Float64}
    lower::Vector{Float64}         # bounds on the variable cells
    upper::Vector{Float64}
    beta::Float64
    jtj::Vector{Float64}           # diag JᵀJ/σ² at the start, for the LSQR column scaling
    options::Inv3DOptions
end

_full(P::Inv3DProblem, x) = (ρ = copy(P.ρ); ρ[P.var] .= x; ρ)
_embed(P::Inv3DProblem, v) = (u = zeros(length(P.ρ)); u[P.var] .= v; u)

predicted(P::Inv3DProblem, x) = forward(P.op, _full(P, x))
residual(P::Inv3DProblem, x) = (predicted(P, x) .- P.d) ./ P.σ
phi_d(P::Inv3DProblem, x) = 0.5sum(abs2, residual(P, x))
rms(P::Inv3DProblem, x) = sqrt(2phi_d(P, x) / length(P.d))
phi_m(P::Inv3DProblem, x) = 0.5P.alpha_s * sum(abs2, sqrt.(P.w) .* (x .- P.xref)) +
                            sum(0.5P.alpha[i] * sum(abs2, P.W[i] * x) for i in 1:3)
objective(P::Inv3DProblem, x) = phi_d(P, x) + P.beta * phi_m(P, x)

reg_hessian(P::Inv3DProblem, v) = P.alpha_s .* P.w .* v .+ sum(P.alpha[i] .* (P.W[i]' * (P.W[i] * v)) for i in 1:3)
data_hessian(P::Inv3DProblem, x, v) = (ρ = _full(P, x); jtvec(P.op, ρ, jvec(P.op, ρ, _embed(P, v)) ./ P.σ .^ 2)[P.var])

function gradient(P::Inv3DProblem, x)
    ρ = _full(P, x)
    g = jtvec(P.op, ρ, (forward(P.op, ρ) .- P.d) ./ P.σ .^ 2)[P.var]
    g .+ P.beta .* (P.alpha_s .* P.w .* (x .- P.xref) .+ sum(P.alpha[i] .* (P.W[i]' * (P.W[i] * x)) for i in 1:3))
end

# first differences between neighbouring variable cells along `dir`, scaled by min(h)/h and √(mean w)
function _differences(m::Model, cells::Vector{CartesianIndex{3}}, var::Vector{Int}, w::Vector{Float64}, dir::Int)
    idx = zeros(Int, size(m.A))
    for (a, c) in enumerate(var)
        idx[cells[c]] = a
    end
    h = (m.dx, m.dy, m.dz)[dir]
    step = CartesianIndex(ntuple(i -> i == dir ? 1 : 0, 3))
    I, J, V, hf, wf = Int[], Int[], Float64[], Float64[], Float64[]
    r = 0
    for c in CartesianIndices(idx)
        c2 = c + step
        checkbounds(Bool, idx, c2) || continue
        a, b = idx[c], idx[c2]
        a > 0 && b > 0 || continue
        r += 1
        push!(I, r, r); push!(J, a, b); push!(V, -1.0, 1.0)
        push!(hf, (h[c[dir]] + h[c2[dir]]) / 2); push!(wf, (w[a] + w[b]) / 2)
    end
    r == 0 && return spzeros(0, length(var))
    Diagonal(sqrt.(wf) .* minimum(hf) ./ hf) * sparse(I, J, V, r, length(var))
end

# diag JᵀJ/σ² = Σ_n G[n, j]² / σ_n² over the stations, for every column j
function _jtj_diag(op::GravForward, σ::Vector{Float64})
    s = 1 ./ σ .^ 2
    if op.solver == :dense
        G = op.G
        out = zeros(size(G, 2))
        Threads.@threads for j in axes(G, 2)
            t = 0.0
            @inbounds for n in axes(G, 1)
                t += s[n] * G[n, j]^2
            end
            out[j] = t
        end
        return out
    end
    _accumulate_rows((a, n, row) -> (a .+= s[n] .* row .^ 2), op)
end

function _weights(op::GravForward, jtj, var, opt::Inv3DOptions)
    opt.weighting == :none && return ones(length(var))
    if opt.weighting == :sensitivity
        w = sqrt.(jtj[var])
    else
        # depth of the cell centre below the mean station height, at least half its thickness
        zs = sum(op.zs) / ndata(op.d)
        w = [max(op.m.cz[op.cells[c][3]] - zs, op.m.dz[op.cells[c][3]] / 2)^(-opt.depth_exponent) for c in var]
    end
    w ./= maximum(w)
    max.(w, 1e-12)
end

"""
    Inv3DProblem(op, options; start=zeros, variable=trues, reference=zeros, lower=bounds[1], upper=bounds[2])

Inverse problem on the forward operator `op` (data and errors from `op.d`). `start`, `reference`,
`lower` and `upper` are active-cell vectors; `variable` marks the inverted active cells. The per-cell
bounds default to `options.bounds`.
"""
function Inv3DProblem(op::GravForward, opt::Inv3DOptions = Inv3DOptions(); start::AbstractVector = zeros(ncells(op)),
                      variable::AbstractVector{Bool} = trues(ncells(op)), reference::AbstractVector = zeros(ncells(op)),
                      lower::AbstractVector = fill(opt.bounds[1], ncells(op)),
                      upper::AbstractVector = fill(opt.bounds[2], ncells(op)))
    _inv3d_validate(opt)
    length(start) == length(variable) == length(reference) == length(lower) == length(upper) == ncells(op) ||
        throw(ArgumentError("start, variable, reference and bounds need one value per active cell"))
    var = findall(variable)
    isempty(var) && throw(ArgumentError("no variable cells"))
    lo, hi = Float64.(lower[var]), Float64.(upper[var])
    all(isfinite, lo) && all(isfinite, hi) && all(lo .<= hi) ||
        throw(ArgumentError("per-cell bounds must be finite with lower ≤ upper"))
    ρ = Float64.(start); ρ[var] .= clamp.(ρ[var], lo, hi)
    jtj = _jtj_diag(op, op.d.error)
    w = _weights(op, jtj, var, opt)
    W = [_differences(op.m, op.cells, var, w, i) for i in 1:3]
    Inv3DProblem(op, copy(op.d.value), copy(op.d.error), var, ρ, Float64.(reference[var]), w, W,
                 [vec(sum(abs2, Wi; dims = 1)) for Wi in W], opt.smallness,
                 (opt.smooth_x, opt.smooth_y, opt.smooth_z), lo, hi, 1.0, jtj[var], opt)
end

# β₀ = ratio × median of vᵀ H_d v / vᵀ H_m v over random v
function _initial_beta(P::Inv3DProblem, x, ratio; n = 5)
    rng = Xoshiro(20261017)
    ratio * median([(v = rand(rng, length(x)); dot(v, data_hessian(P, x, v)) / dot(v, reg_hessian(P, v))) for _ in 1:n])
end

#---------- line search ----------

# armijo backtracking along a capped direction, projected on the bounds; nothing when no step is accepted
function _inv3d_linesearch(P::Inv3DProblem, x, f, direction, gradient)
    opt = P.options
    d = direction .* min(1.0, opt.max_step / max(norm(direction, Inf), eps()))
    for ls in 0:opt.max_linesearch-1
        alpha = 0.5^ls
        candidate = clamp.(x .+ alpha .* d, P.lower, P.upper)
        delta = candidate .- x
        slope = dot(gradient, delta)
        slope < 0 || continue
        trial = objective(P, candidate)
        if isfinite(trial) && trial <= f + 1e-4 * slope
            return candidate, norm(delta, Inf), alpha
        end
    end
    nothing
end

# one step from x at the current β: direction, line search, retries the algorithm asks for
function _inv3d_step(algorithm, work, P::Inv3DProblem, x)
    gradient = inv3d_prepare!(algorithm, work, P, x)
    f = objective(P, x)
    while true
        direction = inv3d_direction(algorithm, work, P, x, gradient)
        found = _inv3d_linesearch(P, x, f, direction, gradient)
        found === nothing || return found
        inv3d_reject!(algorithm, work) || return nothing
    end
end

#---------- driver ----------

"""
    Invert3D(m, d; algorithm=LSQR3DConfig(), options=Inv3DOptions(), fwd=FwdCtrl3D(),
             mask=trues, reference=nothing, lower=nothing, upper=nothing, logfile=nothing)

Invert the data `d` for the density contrast of the variable cells of `m`, minimizing
`φ_d + β φ_m` (see `Inv3D.jl`) within the density bounds.

- `algorithm`: `LSQR3DConfig` or any `AbstractInversion3D`
- `options`: `Inv3DOptions`, shared by every algorithm
- `fwd`: `FwdCtrl3D`, the forward solver
- `mask`: model-shaped Boolean array of the cells that may change; air is always fixed
- `reference`: `Model` the smallness term pulls toward; zero when not given
- `lower`, `upper`: `Model`s of per-cell density bounds (e.g. from a drillhole log), replacing
  `options.bounds` where given
- `logfile`: path of an iteration log, a row appended after every accepted iteration (the
  six-file method writes `0lsqr.log` with a header)

Returns `Inv3DResult`; no other files are written by this method, the six-file method below
writes a run directory.
"""
function Invert3D(m::Model, d::Data; algorithm::AbstractInversion3D = LSQR3DConfig(),
                  options::Inv3DOptions = Inv3DOptions(), fwd::FwdCtrl3D = FwdCtrl3D(),
                  mask::AbstractArray{Bool,3} = trues(size(m.A)), reference::Union{Nothing,Model} = nothing,
                  lower::Union{Nothing,Model} = nothing, upper::Union{Nothing,Model} = nothing,
                  logfile::Union{Nothing,AbstractString} = nothing)
    opt = options
    _inv3d_validate(opt)
    inv3d_validate(algorithm)
    _validate_ctrl(fwd)
    size(mask) == size(m.A) || throw(DimensionMismatch("mask must match the model"))

    #---------- setup ----------
    op = GravForward(m, d; solver = fwd.solver)
    ref = reference === nothing ? zeros(ncells(op)) : model_vector(op, reference)
    lo = lower === nothing ? fill(opt.bounds[1], ncells(op)) : model_vector(op, lower)
    hi = upper === nothing ? fill(opt.bounds[2], ncells(op)) : model_vector(op, upper)
    P = Inv3DProblem(op, opt; start = model_vector(op, m), variable = mask[op.cells], reference = ref, lower = lo, upper = hi)
    x = P.ρ[P.var]
    P.beta = _initial_beta(P, x, opt.beta)
    work = inv3d_init(algorithm, P, x)
    target = opt.target_rms

    #---------- history ----------
    history = NamedTuple[]
    function record(iteration, x, step, alpha)
        r = residual(P, x)
        chi2 = sum(abs2, r)
        entry = (; iteration, objective = chi2 / 2 + P.beta * phi_m(P, x), chi2, rms = sqrt(chi2 / length(r)),
                   regularization = phi_m(P, x), lambda = P.beta, step, alpha, inv3d_info(algorithm, work)...)
        push!(history, entry)
        logfile === nothing || _log_append(logfile, entry)
        opt.verbose && @printf("%s %3d  RMS %.5g  objective %.6g  step %.3g  lambda %.3g\n",
                               uppercase(inv3d_tag(algorithm)), iteration, entry.rms, entry.objective, step, P.beta)
        opt.verbose && flush(stdout)      # batch logs are block buffered otherwise
    end
    record(0, x, 0.0, 0.0)

    #---------- iterations ----------
    reason = :max_iter
    for iteration in 1:opt.max_iter
        if target > 0 && history[end].rms <= target
            reason = :target_rms
            break
        end
        xprev = x
        found = _inv3d_step(algorithm, work, P, xprev)
        if found === nothing
            reason = :line_search_failed
            break
        end
        x, step, alpha = found
        r = rms(P, x)
        if r <= target && r < target * (1 - opt.target_tolerance)
            # overshot: largest β whose step from the previous model still reaches the target. The
            # previous step need not have converged at its β, so β is first raised by the cooling
            # factor until the step stops at or above the target, then log-bisected
            lo, hi = log(P.beta), log(P.beta * opt.cooling)
            within = false
            for _ in 1:6
                P.beta = exp(hi)
                t = _inv3d_step(algorithm, work, P, xprev)
                rt = t === nothing ? Inf : rms(P, t[1])
                rt > target && break
                lo, (x, step, alpha), r = hi, t, rt
                (within = r >= target * (1 - opt.target_tolerance)) && break
                hi += log(opt.cooling)
            end
            for _ in 1:(within ? 0 : 4)
                P.beta = exp((lo + hi) / 2)
                t = _inv3d_step(algorithm, work, P, xprev)
                rt = t === nothing ? Inf : rms(P, t[1])
                if rt <= target
                    lo, (x, step, alpha), r = log(P.beta), t, rt
                    r >= target * (1 - opt.target_tolerance) && break
                else
                    hi = log(P.beta)
                end
            end
            P.beta = exp(lo)
        end
        record(iteration, x, step, alpha)
        inv3d_accept!(algorithm, work, P, xprev, x)
        if r <= target
            reason = :target_rms
            break
        end
        if iteration > 1 && abs(history[end].rms - history[end-1].rms) < 1e-6 * history[end].rms
            reason = :no_progress
            break
        end
        iteration % opt.cooling_interval == 0 && (P.beta /= opt.cooling)
    end

    ρ = _full(P, x)
    model = Model(m.dx, m.dy, m.dz, model_array(op, ρ); origin = m.origin, rotation = m.rotation, datum = m.datum)
    Inv3DResult(model, forward(op, ρ), history[end].rms, history, op.cells[P.var],
                reason ∉ (:max_iter, :line_search_failed), reason, algorithm, P)
end

#---------- file workflow ----------

# shared options and the algorithm config from the inversion control
function _inv3d_from_ctrl(c::InvCtrl3D)
    options = Inv3DOptions(max_iter = c.max_iter, beta = c.lambda, cooling = c.lambda_cooling,
                           cooling_interval = c.cooling_interval, target_rms = c.target_rms,
                           target_tolerance = c.target_tolerance, smallness = c.smallness, smooth_x = c.smooth_x,
                           smooth_y = c.smooth_y, smooth_z = c.smooth_z, weighting = c.weighting,
                           depth_exponent = c.depth_exponent, bounds = c.bounds, max_step = c.max_step,
                           max_linesearch = c.max_linesearch)
    algorithm = LSQR3DConfig(max_iter = c.lsqr_max_iter, tolerance = c.lsqr_tolerance)
    options, algorithm
end

# run_YYYYmmdd_HHMMSS next to the data, suffixed _2, _3 when started in the same second
function _inv3d_run_dir(data_path::AbstractString)
    base = joinpath(dirname(abspath(data_path)), "run_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    dir, k = base, 1
    while ispath(dir)
        k += 1
        dir = "$(base)_$k"
    end
    mkpath(dir)
    dir
end

function _inv3d_open_run(run_dir, data_path, inputs)
    dir = run_dir === nothing ? _inv3d_run_dir(data_path) : mkpath(String(run_dir))
    mkpath(joinpath(dir, "inputs"))
    foreach(p -> cp(p, joinpath(dir, "inputs", basename(p)); force = true), inputs)
    println("Run directory: ", relpath(dir))
    dir
end


# MT-style run log (as 0vfsa3DMT.log): the inputs and controls as "# key = value" lines, one fixed-width row per
# iteration, appended as the run goes so a running inversion can be followed, and the result as "# key = value" lines
_log_stamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
_log_fields(c) = (; (f => getfield(c, f) for f in fieldnames(typeof(c)))...)
function _log_pairs(io, nt::NamedTuple)
    w = maximum(k -> length(String(k)), keys(nt))
    foreach(k -> println(io, "# ", rpad(String(k), w), " = ", nt[k]), keys(nt))
end

function _log_header(path, title, nt::NamedTuple)
    open(path, "w") do io
        println(io, "# ", title, " — ", _log_stamp())
        _log_pairs(io, nt)
    end
    path
end

function _log_footer(path, nt::NamedTuple)
    open(path, "a") do io
        println(io, "# result — ", _log_stamp())
        _log_pairs(io, nt)
    end
    path
end

const _LOG_COLUMN = (rms = "RMS", chi2 = "Chi2", regularization = "PhiM", lsqr = "LSQR")
_log_column(k::Symbol) = get(_LOG_COLUMN, k, titlecase(String(k)))

# one row per iteration, the first field the iteration; the column names before iteration zero
function _log_append(path, entry::NamedTuple)
    cell(v) = v isa Integer ? @sprintf("%14d", v) : @sprintf("%14.6g", v)
    open(path, "a") do io
        if first(entry) == 0
            rule = repeat("-", 6 + 14 * (length(entry) - 1))
            println(io, rule)
            println(io, lpad("Iter", 6), join(lpad.(_log_column.(Base.tail(keys(entry))), 14)))
            println(io, rule)
        end
        println(io, lpad(first(entry), 6), join(cell.(Base.tail(values(entry)))))
    end
end

# predicted data on the observed stations, with the observed errors
_inv3d_predicted(observed::Data, value::AbstractVector) = (p = deepcopy(observed); p.value .= value; p)

# data and cells of the run, in the result of every log
_inv3d_data_label(d::Data) = string(ndata(d), " ", _COMPONENT_NAME[d.component], " (", units(d.component), ")")
_inv3d_cells(m::Model, variable::AbstractArray{Bool,3}) =
    (; mesh = join(size(m.A), " x "), active_cells = count(active_cells(m)), variable_cells = count(variable),
       air_cells = count(air_cells(m)))

"""
    Invert3D(start_path, data_path, fwd_path, inv_path, mask_path, prior_path; run_dir=nothing, bounds=nothing)

Deterministic inversion (LSQR) from six files: the start model, the observed data
(native file, in the mesh frame, see `rotate_data`), `FwdCtrl`, the inversion control, the
mask (0 fixed, 1 variable) and the prior model. The algorithm comes from the inversion
control, the forward solver from `FwdCtrl`. The prior is the reference model of the
smallness term. `bounds = (lower_path, upper_path)` gives per-cell bounds models on the same mesh,
replacing the control's bounds (e.g. a drillhole log).

Everything is written to `run_dir`, by default `run_YYYYmmdd_HHMMSS/` next to the data:
`model.den` (restartable), `data.pred`, the run log `0lsqr.log` ( inputs and controls,
one row per iteration, result) and copies of the inputs in `inputs/`. Returns a named tuple with the run directory, models, data, history, rms and
termination reason.
"""
function Invert3D(start_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString,
                  inv_path::AbstractString, mask_path::AbstractString, prior_path::AbstractString;
                  run_dir::Union{Nothing, AbstractString} = nothing,
                  bounds::Union{Nothing, Tuple{AbstractString, AbstractString}} = nothing)
    fwd, ctrl = ReadFwdCtrl3D(fwd_path), ReadInvCtrl3D(inv_path)
    observed = load_data_grav(data_path)
    start, prior = load_model(start_path), load_model(prior_path)
    same_mesh(start, prior) || error("the prior model must be on the start model's mesh")
    mask = load_mask(mask_path, start)
    any(mask) || error("the mask fixes every cell")
    lower, upper = bounds === nothing ? (nothing, nothing) : load_model.(bounds)
    bounds === nothing || (same_mesh(start, lower) && same_mesh(start, upper)) ||
        error("the bounds models must be on the start model's mesh")
    files = (start_path, data_path, fwd_path, inv_path, mask_path, prior_path)
    dir = _inv3d_open_run(run_dir, data_path, bounds === nothing ? files : (files..., bounds...))

    options, algorithm = _inv3d_from_ctrl(ctrl)
    label = string(ctrl.algorithm)                                                  # run label: log name and titles
    logfile = _log_header(joinpath(dir, "0$label.log"), "$(uppercase(label)) 3D gravity inversion",
                          (; start_model = start_path, data = data_path, fwd_ctrl = fwd_path, inv_ctrl = inv_path,
                             mask = mask_path, prior = prior_path,
                             bounds_models = bounds === nothing ? "none" : join(bounds, " "),
                             solver = fwd.solver, _log_fields(ctrl)...))
    elapsed = @elapsed result = Invert3D(start, observed; algorithm, options, fwd, mask, reference = prior, lower, upper,
                                         logfile)
    predicted = _inv3d_predicted(observed, result.predicted)

    write_model(joinpath(dir, "model.den"), result.model; comment = "$(uppercase(label)) density contrast (kg/m3)")
    write_data_grav(joinpath(dir, "data.pred"), predicted; comment = "predicted")
    _log_footer(logfile, (; termination = result.reason, converged = result.converged, rms = @sprintf("%.6f", result.rms),
                         data = _inv3d_data_label(observed), accepted_iterations = length(result.history) - 1,
                         run_time_s = @sprintf("%.1f", elapsed),
                         final_lambda = @sprintf("%.6g", result.history[end].lambda),
                         density_range = @sprintf("%.6g %.6g", extrema(filter(!isnan, result.model.A))...),
                         _inv3d_cells(start, mask)...))
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, fwd, observed, predicted, start, prior,
       final = result.model, active = mask, history = result.history, rms = result.rms,
       reason = result.reason, converged = result.converged)
end
