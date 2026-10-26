# 3D gravity deterministic inversion framework
# Author: @pankajkmishra
# Algorithm-agnostic problem setup, objective, weighted regularization, projected line search, beta cooling and driver loop
# Each algorithm lives in its own Inv3D_<Name>.jl as a config struct `A <: AbstractInversion3D` plus the
# methods inv3d_init (mutable scratch), inv3d_prepare! (the gradient, once per iteration), inv3d_direction
# (the search direction in density contrast), inv3d_reject! (true = retry with a new direction),
# inv3d_accept!, inv3d_info (a NamedTuple for the history), inv3d_tag (a short log label) and
# inv3d_validate. The driver owns the stopping tests, beta cooling, the bounds and the projected Armijo
# line search, so a new algorithm only decides which direction to take

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using Statistics

#---------- algorithm interface ----------

"""
    AbstractInversion3D

Supertype of 3D deterministic inversion algorithms.
"""
abstract type AbstractInversion3D end

inv3d_reject!(::AbstractInversion3D, work) = false
inv3d_accept!(::AbstractInversion3D, work, problem, old, new) = nothing
inv3d_info(::AbstractInversion3D, work) = (;)
inv3d_validate(::AbstractInversion3D) = nothing

#---------- options and result ----------

"""
    Inv3DOptions(; kwargs...)

Algorithm-independent inversion controls, see `InvCtrl3D` for their meaning.
"""
Base.@kwdef struct Inv3DOptions
    max_iter::Int = 20
    target_rms::Float64 = 1.0
    beta::Union{Nothing, Float64} = nothing
    beta_cooling::Float64 = 2.0
    bounds::Tuple{Float64, Float64} = (-1000.0, 1000.0)
    smallness::Float64 = 0.01
    smooth_x::Float64 = 1.0
    smooth_y::Float64 = 1.0
    smooth_z::Float64 = 1.0
    weighting::Symbol = :sensitivity
    depth_exponent::Float64 = 2.0
    depth_z0::Union{Nothing, Float64} = nothing
    max_linesearch::Int = 12
    step_tolerance::Float64 = 1e-7
    gradient_tolerance::Float64 = 1e-10
    verbose::Bool = true
end

"""
    Inv3DResult

Result of `Invert3D`.
"""
struct Inv3DResult{A<:AbstractInversion3D}
    density::Array{Float64, 3}
    predicted::Vector{Float64}
    chi2::Float64
    rms::Float64
    beta::Float64
    weights::Vector{Float64}
    history::Vector{NamedTuple}
    active_cells::Vector{Int}
    converged::Bool
    reason::Symbol
    algorithm::A
end

#---------- problem ----------

"""
    Inv3DProblem

Fixed part of an inversion: the transposed sensitivity of the active cells, the data and their errors, the weighted regularization operator, the reference model and the options.
"""
struct Inv3DProblem{M<:AbstractMatrix{Float32}}
    Gt::M                                    # (n_active x n_data), mGal per kg/m3
    sigma::Vector{Float64}                   # data errors, mGal
    dw::Vector{Float64}                      # (observed - fixed cells) / sigma
    cells::Vector{Int}                       # active cells, linear model indices
    R::SparseMatrixCSC{Float64, Int}         # weighted regularizer on the active cells
    RtR::SparseMatrixCSC{Float64, Int}
    mref::Vector{Float64}
    weights::Vector{Float64}
    options::Inv3DOptions
end

function _inv3d_validate(opt::Inv3DOptions)
    opt.max_iter >= 0 || throw(ArgumentError("max_iter must be nonnegative"))
    opt.beta === nothing || (isfinite(opt.beta) && opt.beta >= 0) || throw(ArgumentError("beta must be nothing or nonnegative"))
    opt.beta_cooling >= 1 || throw(ArgumentError("beta_cooling must be at least 1"))
    opt.bounds[1] < opt.bounds[2] || throw(ArgumentError("bounds must be increasing"))
    all(x -> isfinite(x) && x >= 0, (opt.smallness, opt.smooth_x, opt.smooth_y, opt.smooth_z, opt.target_rms,
                                     opt.depth_exponent, opt.step_tolerance, opt.gradient_tolerance)) ||
        throw(ArgumentError("regularization and stopping controls must be finite and nonnegative"))
    opt.weighting in (:sensitivity, :depth, :none) || throw(ArgumentError("weighting must be :sensitivity, :depth or :none"))
    opt.max_linesearch > 0 || throw(ArgumentError("max_linesearch must be positive"))
    nothing
end

# J v and Jt w of the data-weighted sensitivity, J = diag(1/sigma) G; the matrix products run in Float32
_inv3d_J(p::Inv3DProblem, v) = Float64.(p.Gt' * Float32.(v)) ./ p.sigma
_inv3d_Jt(p::Inv3DProblem, w) = Float64.(p.Gt * Float32.(w ./ p.sigma))

# diag(JtJ), the squared column norms of the weighted sensitivity
function _inv3d_JtJ_diag(Gt, sigma)
    w2  = Float32.(1 ./ sigma .^ 2)
    out = zeros(size(Gt, 1))
    Threads.@threads for j in axes(Gt, 1)
        s = 0.0
        @inbounds for i in axes(Gt, 2)
            s += Float64(Gt[j, i])^2 * w2[i]
        end
        out[j] = s
    end
    out
end

#---------- model weighting ----------

# depth of each active cell centre below the top of the earth in its column
function _inv3d_cell_depths(m::GravModel3D, cells)
    surface = surface_depth(m)
    nx, ny, _ = size(m.A)
    [begin
        i, j, k = Tuple(CartesianIndices(m.A)[c])
        m.cz[k] - surface[i, j]
     end for c in cells]
end

# z0 of the depth weighting: median sensor height above the ground plus half the top cell
function _inv3d_auto_z0(m::GravModel3D, d::GravData3D)
    surface = surface_depth(m)
    nx, ny, _ = size(m.A)
    top   = minimum(filter(isfinite, surface))
    depth = sensor_depth(m, d)
    h     = map(eachindex(d.x)) do i
        ix = clamp(searchsortedlast(m.x, d.x[i]), 1, nx); iy = clamp(searchsortedlast(m.y, d.y[i]), 1, ny)
        s  = surface[ix, iy]
        (isfinite(s) ? s : top) - depth[i]
    end
    max(median(h), 0.0) + m.dz[1] / 2
end

"""
    model_weights(model, data, cells, opt[, colnorm2]) -> Vector

Model weights of the active cells, normalised to a maximum of 1.
"""
function model_weights(m::GravModel3D, d::GravData3D, cells, opt::Inv3DOptions,
                       colnorm2::Union{Nothing, AbstractVector} = nothing)
    opt.weighting == :none && return ones(length(cells))
    w = if opt.weighting == :depth
        z0 = opt.depth_z0 === nothing ? _inv3d_auto_z0(m, d) : opt.depth_z0
        (max.(_inv3d_cell_depths(m, cells), 0.0) .+ z0) .^ (-opt.depth_exponent / 2)
    else
        colnorm2 === nothing && throw(ArgumentError("sensitivity weighting needs the column norms"))
        vol = vec(cell_volumes(m))[cells]
        sqrt.(sqrt.(colnorm2) ./ (vol ./ median(vol)))
    end
    w ./= maximum(w)
    max.(w, 1e-6)
end

#---------- regularizer ----------

# unweighted finite-volume regularizer of the active cells; fixed and air cells take no term and no pair
function _inv3d_regularizer(m::GravModel3D, cells, opt::Inv3DOptions)
    nx, ny, nz = size(m.A)
    index = zeros(Int, nx, ny, nz)
    for (n, c) in enumerate(cells)
        index[c] = n
    end
    V0 = median(m.dx) * median(m.dy) * median(m.dz)
    h0 = cbrt(V0)
    I, J, V = Int[], Int[], Float64[]
    row = 0
    pair!(a, b, w) = (row += 1; append!(I, (row, row)); append!(J, (a, b)); append!(V, (-w, w)))
    for k in 1:nz, j in 1:ny, i in 1:nx
        a = index[i, j, k]
        a == 0 && continue
        if opt.smallness > 0
            row += 1
            push!(I, row); push!(J, a); push!(V, sqrt(opt.smallness * m.dx[i] * m.dy[j] * m.dz[k] / V0))
        end
        if i < nx && opt.smooth_x > 0 && (b = index[i+1, j, k]) != 0
            pair!(a, b, sqrt(opt.smooth_x * m.dy[j] * m.dz[k] / ((m.dx[i] + m.dx[i+1]) / 2) / h0))
        end
        if j < ny && opt.smooth_y > 0 && (b = index[i, j+1, k]) != 0
            pair!(a, b, sqrt(opt.smooth_y * m.dx[i] * m.dz[k] / ((m.dy[j] + m.dy[j+1]) / 2) / h0))
        end
        if k < nz && opt.smooth_z > 0 && (b = index[i, j, k+1]) != 0
            pair!(a, b, sqrt(opt.smooth_z * m.dx[i] * m.dy[j] / ((m.dz[k] + m.dz[k+1]) / 2) / h0))
        end
    end
    sparse(I, J, V, row, length(cells))
end

#---------- state: one evaluated model ----------

function _inv3d_evaluate(problem::Inv3DProblem, m::AbstractVector, beta::Real)
    r   = _inv3d_J(problem, m) .- problem.dw
    reg = problem.R * (m .- problem.mref)
    phid, phim = sum(abs2, r), sum(abs2, reg)
    (; m = Vector{Float64}(m), r, reg, beta = Float64(beta), phid, phim, objective = (phid + beta * phim) / 2)
end

_inv3d_rebeta(problem, state, beta) = (; state..., beta = Float64(beta), objective = (state.phid + beta * state.phim) / 2)

#---------- derivatives shared by all algorithms ----------

"""
    inv3d_gradient(problem, state)

Gradient of the objective in active density contrast, `Jᵀr + βRᵀR(m - mref)`.
"""
inv3d_gradient(problem::Inv3DProblem, state) =
    _inv3d_Jt(problem, state.r) .+ state.beta .* (problem.R' * state.reg)

"""
    inv3d_hessian_mul(problem, state, v)

Gauss–Newton Hessian times a vector, `JᵀJv + βRᵀRv`, without forming the matrix.
"""
inv3d_hessian_mul(problem::Inv3DProblem, state, v) =
    _inv3d_Jt(problem, _inv3d_J(problem, v)) .+ state.beta .* (problem.RtR * v)

# cells held at a bound, where the negative gradient points out of the box
function inv3d_bound_cells(problem::Inv3DProblem, state, gradient)
    lo, hi = problem.options.bounds
    tol = 1e-10 * max(1.0, hi - lo)
    [(m <= lo + tol && g > 0) || (m >= hi - tol && g < 0) for (m, g) in zip(state.m, gradient)]
end

# largest eigenvalue of a symmetric operator by power iteration
function _inv3d_power(mul, n; iterations = 20)
    v      = normalize!(ones(n) .+ 0.1 .* sin.(1:n))
    lambda = 0.0
    for _ in 1:iterations
        w      = mul(v)
        lambda = norm(w)
        lambda == 0 && return 0.0
        v = w ./ lambda
    end
    lambda
end

#---------- line search ----------

# projected armijo backtracking: the trial is clipped to the bounds; nothing when no step is accepted
function _inv3d_linesearch(problem::Inv3DProblem, state, direction, gradient)
    opt = problem.options
    lo, hi = opt.bounds
    for ls in 0:opt.max_linesearch-1
        alpha     = 0.5^ls
        candidate = clamp.(state.m .+ alpha .* direction, lo, hi)
        delta     = candidate .- state.m
        slope     = dot(gradient, delta)
        slope < 0 || continue
        trial = _inv3d_evaluate(problem, candidate, state.beta)
        if isfinite(trial.objective) && trial.objective <= state.objective + 1e-4 * slope
            return trial, norm(delta, Inf), alpha
        end
    end
    nothing
end

#---------- target misfit ----------

# the cooled step overshot the target: bisect beta geometrically for an rms in [0.9, 1] x target
function _inv3d_refine_beta(algorithm, work, problem, old, found, target, nd; attempts::Int = 4)
    best = found
    lo, hi = old.beta, old.beta * problem.options.beta_cooling
    for _ in 1:attempts
        beta      = sqrt(lo * hi)
        trial_old = _inv3d_rebeta(problem, old, beta)
        gradient  = inv3d_prepare!(algorithm, work, problem, trial_old)
        direction = inv3d_direction(algorithm, work, problem, trial_old, gradient)
        trial     = _inv3d_linesearch(problem, trial_old, direction, gradient)
        trial === nothing && break
        r = sqrt(trial[1].phid / nd)
        if r <= target
            best = trial
            r >= 0.9 * target && break
            lo = beta
        else
            hi = beta
        end
    end
    best
end

#---------- driver ----------

"""
    Invert3D(model, observed::GravData3D, fwd::FwdCtrl3D; algorithm=GaussNewton3DConfig(),
             options=Inv3DOptions(), active_cells=nothing, reference=nothing) -> Inv3DResult

Invert gravity data for density contrast on the model's mesh, minimizing
"""
function Invert3D(model::GravModel3D, observed::GravData3D, fwd::FwdCtrl3D;
                  algorithm::AbstractInversion3D = GaussNewton3DConfig(),
                  options::Inv3DOptions = Inv3DOptions(),
                  active_cells = nothing, reference::Union{Nothing, GravModel3D} = nothing)
    opt = options
    _inv3d_validate(opt)
    inv3d_validate(algorithm)
    lo, hi = opt.bounds

    #---------- setup ----------
    air   = air_mask(model)
    cells = if active_cells === nothing
        findall(vec(.!air))
    elseif active_cells isa AbstractArray{Bool}
        size(active_cells) == size(model.A) || throw(DimensionMismatch("active mask must match the model"))
        findall(vec(active_cells))
    else
        collect(Int, active_cells)
    end
    isempty(cells) && throw(ArgumentError("at least one active earth cell is required"))
    any(c -> air[c], cells) && throw(ArgumentError("air cells cannot be inverted"))
    ref = reference === nothing ? model : reference
    same_mesh(ref, model) || throw(ArgumentError("the reference model must be on the start model's mesh"))
    sigma = copy(observed.error)
    all(s -> isfinite(s) && s > 0, sigma) || throw(ArgumentError("data errors must be positive and finite"))

    # fixed earth cells keep their start values; their field is taken off the data once
    fixed = falses(length(model.A))
    fixed[findall(vec(.!air))] .= true
    fixed[cells] .= false
    fixed_model = model_like(model, ifelse.(reshape(fixed, size(model.A)), model.A, 0.0))
    dfixed      = any(fixed) && any(!=(0), filter(!isnan, fixed_model.A)) ? forward3d(fixed_model, observed, fwd) : zeros(ndata(observed))

    opt.verbose && @printf("Inversion: %d data, %d active of %d cells, %d threads\n",
                           ndata(observed), length(cells), length(model.A), Threads.nthreads())
    tG = @elapsed Gt = Sensitivity3D(model, observed, fwd; cells, verbose = opt.verbose)
    opt.verbose && @printf("Sensitivity: %.1f s\n", tG)
    colnorm2 = _inv3d_JtJ_diag(Gt, sigma)
    weights  = model_weights(model, observed, cells, opt, colnorm2)
    R        = _inv3d_regularizer(model, cells, opt) * Diagonal(weights)
    problem = Inv3DProblem(Gt, sigma, (observed.value .- dfixed) ./ sigma, cells, R, sparse(R' * R),
                           clamp.(replace(ref.A[cells], NaN => 0.0), lo, hi), weights, opt)

    beta0 = if opt.beta === nothing
        λd = _inv3d_power(v -> _inv3d_Jt(problem, _inv3d_J(problem, v)), length(cells))
        λm = _inv3d_power(v -> problem.RtR * v, length(cells))
        λd / max(λm, eps())
    else
        opt.beta
    end
    m0    = clamp.(replace(model.A[cells], NaN => 0.0), lo, hi)
    state = _inv3d_evaluate(problem, m0, beta0)
    all(isfinite, state.r) || error("initial forward response is not finite")
    work = inv3d_init(algorithm, problem, state)

    #---------- history ----------
    history = NamedTuple[]
    nd      = ndata(observed)
    rms_of(s) = sqrt(s.phid / nd)
    function record(iteration, step, alpha, info)
        entry = (; iteration, objective = state.objective, chi2 = state.phid, rms = sqrt(state.phid / nd),
                   regularization = state.phim / 2, beta = state.beta, step, alpha,
                   at_bounds = count(x -> x <= lo || x >= hi, state.m), info...)
        push!(history, entry)
        opt.verbose && @printf("%s %3d  RMS %.5g  beta %.4g  objective %.6g  step %.3g\n",
                               uppercase(inv3d_tag(algorithm)), iteration, entry.rms, entry.beta, entry.objective, step)
        opt.verbose && flush(stdout)      # batch logs are block buffered otherwise
    end
    record(0, 0.0, 0.0, inv3d_info(algorithm, work))

    #---------- iterations ----------
    reached() = opt.target_rms > 0 && history[end].rms <= opt.target_rms
    reason = :max_iter
    for iteration in 1:opt.max_iter
        if reached()
            reason = :target_rms
            break
        end
        gradient  = inv3d_prepare!(algorithm, work, problem, state)
        projected = ifelse.(inv3d_bound_cells(problem, state, gradient), 0.0, gradient)
        if norm(projected, Inf) <= opt.gradient_tolerance
            reason = :gradient_tolerance
            break
        end

        # the algorithm may retry with a new direction (e.g. more damping) after a failed search
        found = nothing
        while true
            direction = inv3d_direction(algorithm, work, problem, state, gradient)
            found     = _inv3d_linesearch(problem, state, direction, gradient)
            found === nothing || break
            inv3d_reject!(algorithm, work) || break
        end
        if found === nothing
            reason = :line_search_failed
            break
        end

        old = state
        state, step, alpha = found
        if opt.target_rms > 0 && opt.beta_cooling > 1 && iteration > 1 && rms_of(state) < 0.9 * opt.target_rms
            state, step, alpha = _inv3d_refine_beta(algorithm, work, problem, old, (state, step, alpha), opt.target_rms, nd)
        end
        record(iteration, step, alpha, inv3d_info(algorithm, work))
        inv3d_accept!(algorithm, work, problem, old, state)

        # cool beta while the data are not yet fitted; with a fixed beta a vanishing step ends the run
        if !reached() && opt.beta_cooling > 1
            state = _inv3d_rebeta(problem, state, state.beta / opt.beta_cooling)
        elseif step <= opt.step_tolerance
            reason = :step_tolerance
            break
        end
    end
    reached() && (reason = :target_rms)

    rho        = copy(model.A)
    rho[cells] = state.m
    pred       = _inv3d_J(problem, state.m) .* sigma .+ dfixed
    Inv3DResult(rho, pred, state.phid, sqrt(state.phid / nd), state.beta, weights, history, cells,
                reason ∉ (:max_iter, :line_search_failed), reason, algorithm)
end

#---------- file workflow ----------

# shared options and the algorithm config from InvCtrl
function _inv3d_from_ctrl(c::InvCtrl3D)
    options = Inv3DOptions(max_iter = c.max_iter, target_rms = c.target_rms, beta = c.lambda,
                           beta_cooling = c.lambda_cooling, bounds = c.bounds, smallness = c.smallness,
                           smooth_x = c.smooth_x, smooth_y = c.smooth_y, smooth_z = c.smooth_z,
                           weighting = c.weighting, depth_exponent = c.depth_exponent, depth_z0 = c.depth_z0,
                           max_linesearch = c.max_linesearch)
    algorithm = GaussNewton3DConfig(damping = c.gn_damping, cg_iterations = c.cg_iterations, cg_tolerance = c.cg_tolerance)
    options, algorithm
end

# the error used: the largest of the file error, the floor and the percentage of |datum|
function _inv3d_errors(d::GravData3D, c::Union{InvCtrl3D, InvCtrlINR3D})
    e = max.(d.error, c.error_floor, c.error_percent / 100 .* abs.(d.value))
    all(>(0), e) || error("some data errors are zero; set 'Error floor (mGal)' in the inversion control")
    data_like(d, d.value; error = e)
end

# <label>_YYYYmmdd_HHMMSS, suffixed _2, _3 when two runs start in the same second
function _inv3d_run_dir(data_path::AbstractString, label = nothing)
    stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    base  = label === nothing ? joinpath(dirname(abspath(data_path)), "run_" * stamp) :
                               rstrip(String(label), ('_', '/', '\\')) * "_" * stamp
    dir, k = base, 1
    while ispath(dir)
        k += 1
        dir = "$(base)_$k"
    end
    mkpath(dir)
    dir
end

function _inv3d_open_run(run_dir, data_path, inputs)
    dir = _inv3d_run_dir(data_path, run_dir)
    mkpath(joinpath(dir, "inputs"))
    foreach(p -> cp(p, joinpath(dir, "inputs", basename(p)); force = true), inputs)
    println("Run directory: ", dir)
    dir
end

_inv3d_hist_cell(v::Integer)       = string(v)
_inv3d_hist_cell(v::AbstractFloat) = @sprintf("%.6g", v)
_inv3d_hist_cell(v)                = string(v)

# a plain-text, fixed-width convergence log: a title line, then the column table between
# dashed rules, each column as wide as its header or its widest value
function _inv3d_write_history(path, history, title)
    names = string.(keys(first(history)))
    rows  = [_inv3d_hist_cell.(values(h)) for h in history]
    width(i) = max(length(names[i]), maximum(r -> length(r[i]), rows))
    widths = [width(i) for i in eachindex(names)]
    rule   = "-"^(sum(widths) + length(widths) - 1)
    row(cells) = join((lpad(cells[i], widths[i]) for i in eachindex(cells)), " ")
    open(path, "w") do io
        println(io, "# ", title, " — ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, rule)
        println(io, row(names))
        println(io, rule)
        foreach(r -> println(io, row(r)), rows)
        println(io, rule)
    end
    path
end

"""
    Invert3D(start_path, data_path, fwd_path, inv_path, cov_path, prior_path; run_dir=nothing)

Inversion from six files: the start model, the observed data, `FwdCtrl`, an inversion control file, the model covariance and the prior model.
"""
function Invert3D(start_path::AbstractString, data_path::AbstractString, fwd_path::AbstractString,
                  inv_path::AbstractString, cov_path::AbstractString, prior_path::AbstractString;
                  run_dir::Union{Nothing, AbstractString} = nothing)
    fwd, ctrl, cov = ReadFwdCtrl3D(fwd_path), ReadInvCtrl(inv_path), ReadCov3D(cov_path)
    start, prior = ReadModel3D(start_path), ReadModel3D(prior_path)
    same_mesh(start, prior) || error("the prior model must be on the start model's mesh")
    size(cov.mask) == size(start.A) || error("the covariance mask is $(size(cov.mask)) cells, the model $(size(start.A))")
    air = air_mask(start)
    bad = findall(air .& (cov.mask .!= 0))
    isempty(bad) || error("$(length(bad)) air cell(s) have covariance mask ≠ 0, first at $(Tuple(bad[1]))")
    active = (cov.mask .!= 0) .& .!air
    any(active) || error("the covariance mask fixes every cell")
    if !ctrl.invert_padding
        # the lateral padding is held fixed on top of the mask; every layer stays free
        ix, iy = lateral_core_ranges(start)
        keep = falses(size(active))
        keep[ix, iy, :] .= true
        active .&= keep
        any(active) || error("'Invert padding : no' leaves no free cell: the covariance mask fixes the whole core")
    end
    observed = _inv3d_errors(ReadData3D(data_path), ctrl)

    dir   = _inv3d_open_run(run_dir, data_path, (start_path, data_path, fwd_path, inv_path, cov_path, prior_path))
    is_gn = ctrl isa InvCtrl3D
    if is_gn
        options, algorithm = _inv3d_from_ctrl(ctrl)
        result = Invert3D(start, observed, fwd; algorithm, options, active_cells = active, reference = prior)
    else
        options, config = _inr_from_ctrl(ctrl)
        result = INR3D(start, observed, fwd; config, options, active_cells = active)
    end
    GC.gc()
    if fwd.write_frechet
        Gt = Sensitivity3D(start, observed, fwd; cells = result.active_cells)
        SaveSensitivity3D(joinpath(dir, "sensitivity.jld2"), Gt; cells = result.active_cells)
    end

    tag       = uppercase(string(ctrl.algorithm))
    final     = model_like(start, result.density; name = joinpath(dir, "model.den"))
    predicted = data_like(observed, result.predicted; title = "Predicted $(observed.title)")
    WriteModel3D(joinpath(dir, "model.den"), final; title = "$tag inversion of $(basename(data_path))")
    WriteData3D(joinpath(dir, "data.pred"), predicted)
    _inv3d_write_history(joinpath(dir, "History.log"), result.history, "$tag 3D gravity inversion")
    open(joinpath(dir, "Summary.txt"), "w") do io
        println(io, "Algorithm: ", tag)
        println(io, "Termination: ", result.reason)
        println(io, "Converged: ", result.converged)
        @printf(io, "RMS: %.6f\n", result.rms)
        if is_gn
            @printf(io, "Final lambda: %.6g\n", result.beta)
            println(io, "Model weighting: ", ctrl.weighting)
        else
            println(io, "Model weighting: none (implicit, network parameterization)")
            println(io, "Network: width ", ctrl.width, ", layers ", ctrl.layers,
                    ", frequencies ", ctrl.frequencies, ", seed ", ctrl.seed)
        end
        println(io, "Data count: ", ndata(observed))
        println(io, "Accepted iterations: ", length(result.history) - 1)
        println(io, "Active cells: ", length(result.active_cells))
        println(io, "Padding inverted: ", ctrl.invert_padding)
        println(io, "Air cells: ", count(air))
        println(io, "Fixed earth cells: ", count(.!air) - length(result.active_cells))
        rng = extrema(filter(!isnan, final.A))
        @printf(io, "Density contrast range: %.4g to %.4g kg/m3\n", rng...)
    end
    (; run_dir = dir, algorithm = ctrl.algorithm, ctrl, fwd, observed, predicted, start, prior, final,
       active, weights = result.weights, history = result.history, rms = result.rms, beta = result.beta,
       reason = result.reason, converged = result.converged)
end
