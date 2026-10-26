# 3D gravity implicit-neural-representation (INR) inversion
# Author: @pankajkmishra
# The density contrast of each active cell is the (bounded) output of a coordinate MLP
# (Lux.jl): rho(x) = mid + half * tanh(f_theta(gamma(x))), gamma a per-axis-normalized
# Fourier/positional feature map of the cell centre. Only the data misfit is minimised,
# by Adam, stopped at the target rms (discrepancy principle); there is no explicit
# smoothness or depth-weighting term -- the network's limited capacity plus early
# stopping are the implicit regularization.
#
# Unlike Inv3D_GN.jl, this algorithm does NOT plug into Invert3D's inv3d_direction/line
# search: that driver clamps and backtracks in model space, so an accepted step is no
# longer f_theta(.) for any theta once a step is rejected or clipped to the bounds -- the
# network and the model would drift apart. INR3D() below is a self-contained driver with
# its own loop over theta, sharing only the data term (Inv3DProblem/_inv3d_evaluate/
# inv3d_gradient from Inv3D.jl, with regularization weights forced to zero so the shared
# machinery reduces to the pure data misfit) and the Inv3DResult type, so plotting and the
# six-file run workflow are common to GN and INR.
#
# The gradient of the loss in theta is computed by chaining an *analytic* data-misfit
# gradient in m (Jᵀr, from Inv3D.jl, so Zygote never traces the sensitivity matrix G) with
# a Zygote pullback of the network only (gamma(x) -> f_theta -> rho): d(loss)/d(theta) =
# (d(rho)/d(theta))ᵗ (d(loss)/d(rho)). This is both the efficient and the numerically
# sound choice: G is a large dense/sparse operator, the network is a few hundred cells wide.

using Lux
using Zygote
using Optimisers
using Random

#---------- configuration ----------

"""
    INR3DConfig(; width=256, layers=4, frequencies=2, learning_rate=1e-3, scaling=:per_axis,
                  seed=42, divergence_factor=3.0, max_restarts=3, schedule=:cosine,
                  rate_floor=0.05, anneal=0.5)

Input of `INR3D()`: the network and its training schedule.
"""
Base.@kwdef struct INR3DConfig <: AbstractInversion3D
    width::Int                = 256      # neurons per hidden layer; capacity, and so the regularization
    layers::Int               = 4        # number of hidden layers
    frequencies::Int          = 2        # Fourier feature octaves of the cell centre; more = sharper structure
    learning_rate::Float64    = 1e-3     # peak Adam step size
    scaling::Symbol           = :per_axis # coordinate normalization, :per_axis or :isotropic
    seed::Int                 = 42       # network initialization; change it to test robustness
    divergence_factor::Float64 = 3.0     # restart when the rms exceeds this multiple of the best; 0 disables
    max_restarts::Int         = 3        # how often a diverging run may fall back and halve the rate
    schedule::Symbol          = :cosine  # learning rate decay, :cosine or :none
    rate_floor::Float64       = 0.05     # cosine decays the rate to this fraction of the peak
    anneal::Float64           = 0.5      # fraction of the budget over which Fourier bands are released
end

const INR3DResult = Inv3DResult{INR3DConfig}

inv3d_tag(::INR3DConfig) = "inr"

function inv3d_validate(alg::INR3DConfig)
    alg.width > 0 && alg.layers > 0 && alg.frequencies >= 0 ||
        throw(ArgumentError("hidden width and layers must be positive and the frequencies nonnegative"))
    alg.learning_rate > 0 || throw(ArgumentError("learning rate must be positive"))
    alg.divergence_factor == 0 || alg.divergence_factor > 1 ||
        throw(ArgumentError("divergence factor must be zero or greater than one"))
    alg.max_restarts >= 0 || throw(ArgumentError("the number of restarts must be nonnegative"))
    alg.schedule in (:cosine, :none) || throw(ArgumentError("schedule must be :cosine or :none"))
    0 < alg.rate_floor <= 1 || throw(ArgumentError("the rate floor must be in (0, 1]"))
    0 <= alg.anneal <= 1 || throw(ArgumentError("the annealing fraction must be in [0, 1]"))
    nothing
end

#---------- training schedules ----------

# peak rate after a 2 % warm-up, then a cosine decay to `rate_floor` of the peak at the last epoch
function _inr_rate(alg::INR3DConfig, peak::Real, epoch::Integer, max_iter::Integer)
    alg.schedule === :none && return peak
    warm = max(1, cld(max_iter, 50))
    epoch <= warm && return peak * epoch / warm
    t = (epoch - warm) / max(1, max_iter - warm)
    peak * (alg.rate_floor + (1 - alg.rate_floor) * (1 + cos(pi * t)) / 2)
end

# per-row weights of the feature matrix, releasing the Fourier bands lowest first over `anneal`
function _inr_band_weights(alg::INR3DConfig, nrow::Integer, epoch::Integer, max_iter::Integer)
    w = ones(Float32, nrow)
    (alg.anneal <= 0 || alg.frequencies == 0) && return w
    progress = epoch / max(1, max_iter)
    for k in 0:alg.frequencies-1
        a = clamp(progress * alg.frequencies / alg.anneal - k, 0, 1)
        w[4+6k:min(nrow, 9+6k)] .= Float32(a)         # the sin and cos rows of band k
    end
    w
end

#---------- coordinate features ----------

# cell-centre coordinates of the active cells, standardized per axis (or isotropically) to ~N(0,1)
function _inr_coordinates(m::GravModel3D, cells; scaling::Symbol = :per_axis)
    xs = Vector{Float64}(undef, length(cells)); ys = similar(xs); zs = similar(xs)
    for (n, c) in enumerate(cells)
        i, j, k = Tuple(CartesianIndices(m.A)[c])
        xs[n], ys[n], zs[n] = m.cx[i], m.cy[j], m.cz[k]
    end
    xm, ym, zm = mean(xs), mean(ys), mean(zs)
    sx, sy, sz = if scaling === :per_axis
        max(std(xs), eps()), max(std(ys), eps()), max(std(zs), eps())
    elseif scaling === :isotropic
        s = max(std(xs), std(ys), std(zs), eps())
        s, s, s
    else
        throw(ArgumentError("scaling must be :per_axis or :isotropic"))
    end
    coords = Matrix{Float32}(undef, 3, length(cells))
    coords[1, :] .= Float32.((xs .- xm) ./ sx)
    coords[2, :] .= Float32.((ys .- ym) ./ sy)
    coords[3, :] .= Float32.((zs .- zm) ./ sz)
    coords
end

# Fourier/positional features: [x; sin(x); cos(x); sin(2x); cos(2x); ...]
function _inr_features(coords::AbstractMatrix{Float32}, nfreq::Int)
    feats = Any[coords]
    for k in 0:nfreq-1
        scale = Float32(2.0^k)
        push!(feats, sin.(scale .* coords))
        push!(feats, cos.(scale .* coords))
    end
    vcat(feats...)
end

#---------- network ----------

_inr_small_dense(dim_in, dim_out) =
    Lux.Dense(dim_in => dim_out; init_weight = (rng, dims...) -> 0.01f0 .* Lux.glorot_uniform(rng, dims...))

function _inr_network(input_dim::Int, alg::INR3DConfig)
    hidden = [Lux.Dense(input_dim => alg.width, Lux.leakyrelu)]
    for _ in 2:alg.layers
        push!(hidden, Lux.Dense(alg.width => alg.width, Lux.leakyrelu))
    end
    Lux.Chain(hidden..., _inr_small_dense(alg.width, 1))
end

#---------- entry point ----------

"""
    INR3D(model, observed::GravData3D, fwd::FwdCtrl3D; config=INR3DConfig(),
          options=Inv3DOptions(), active_cells=nothing) -> INR3DResult

Implicit-neural-representation density inversion: the density of each active cell is the bounded output of a coordinate MLP (`options.bounds`), trained by Adam to minimise the data misfit only (`options.target_rms`/`options.max_iter` set the discrepancy-principle early stop).
"""
function INR3D(model::GravModel3D, observed::GravData3D, fwd::FwdCtrl3D;
               config::INR3DConfig = INR3DConfig(), options::Inv3DOptions = Inv3DOptions(),
               active_cells = nothing)
    inv3d_validate(config)
    lo, hi = options.bounds
    mid, half = Float32((lo + hi) / 2), Float32((hi - lo) / 2)

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
    sigma = copy(observed.error)
    all(s -> isfinite(s) && s > 0, sigma) || throw(ArgumentError("data errors must be positive and finite"))

    # fixed earth cells keep their start values; their field is taken off the data once (as in Invert3D/GN)
    fixed = falses(length(model.A))
    fixed[findall(vec(.!air))] .= true
    fixed[cells] .= false
    fixed_model = model_like(model, ifelse.(reshape(fixed, size(model.A)), model.A, 0.0))
    dfixed      = any(fixed) && any(!=(0), filter(!isnan, fixed_model.A)) ? forward3d(fixed_model, observed, fwd) : zeros(ndata(observed))

    options.verbose && @printf("INR: %d data, %d active of %d cells, %d threads\n",
                               ndata(observed), length(cells), length(model.A), Threads.nthreads())
    tG = @elapsed Gt = Sensitivity3D(model, observed, fwd; cells, verbose = options.verbose)
    options.verbose && @printf("Sensitivity: %.1f s\n", tG)

    # zero regularization, so the shared machinery reduces to the pure data misfit and G is never traced
    zero_opt = Inv3DOptions(bounds = options.bounds, smallness = 0.0, smooth_x = 0.0, smooth_y = 0.0,
                            smooth_z = 0.0, weighting = :none, verbose = false)
    R = _inv3d_regularizer(model, cells, zero_opt)
    problem = Inv3DProblem(Gt, sigma, (observed.value .- dfixed) ./ sigma, cells, R, sparse(R' * R),
                           zeros(length(cells)), ones(length(cells)), zero_opt)

    coords   = _inr_coordinates(model, cells; scaling = config.scaling)
    features = _inr_features(coords, config.frequencies)
    rng      = Random.default_rng()
    Random.seed!(rng, config.seed)
    network = _inr_network(size(features, 1), config)
    ps, st = Lux.setup(rng, network)
    opt_state = Optimisers.setup(Optimisers.Adam(config.learning_rate), ps)

    # `X` is the feature matrix of the epoch: the same features, with the Fourier bands faded in
    density_of(p, X) = mid .+ half .* tanh.(vec(first(Lux.apply(network, X, p, st))))
    density_of(p) = density_of(p, features)

    nd    = ndata(observed)
    state = _inv3d_evaluate(problem, Float64.(density_of(ps)), 0.0)
    all(isfinite, state.r) || error("initial network output gives a non-finite forward response")

    history = NamedTuple[]
    push!(history, (; iteration = 0, objective = state.objective, chi2 = state.phid, rms = sqrt(state.phid / nd),
                     regularization = 0.0, beta = 0.0, step = 0.0, alpha = 0.0))
    options.verbose && @printf("INR   0  RMS %.5g\n", history[end].rms)
    flush(stdout)

    # Adam is not monotone, so the best network is kept and the divergence test arms after a warm-up
    best_rms, best_m, best_ps = history[end].rms, Float64.(density_of(ps)), deepcopy(ps)
    warmup = max(25, options.max_iter ÷ 20)
    rate, restarts, last_restart = config.learning_rate, 0, 0
    reason = :max_iter
    for epoch in 1:options.max_iter
        X = features .* _inr_band_weights(config, size(features, 1), epoch, options.max_iter)
        m_theta, back = Zygote.pullback(p -> density_of(p, X), ps)
        state = _inv3d_evaluate(problem, Float64.(m_theta), 0.0)
        rms   = sqrt(state.phid / nd)
        push!(history, (; iteration = epoch, objective = state.objective, chi2 = state.phid, rms,
                         regularization = 0.0, beta = 0.0, step = 0.0, alpha = 0.0))
        if options.verbose && (epoch % max(1, options.max_iter ÷ 10) == 0 || epoch == options.max_iter)
            @printf("INR %3d  RMS %.5g\n", epoch, rms)
            flush(stdout)
        end
        if rms < best_rms
            best_rms, best_m, best_ps = rms, Float64.(m_theta), deepcopy(ps)
        end

        if options.target_rms > 0 && rms <= options.target_rms
            reason = :target_rms
            break
        end
        # a blown-up misfit means too long a step: fall back to the best network and halve the rate
        if config.divergence_factor > 0 && epoch > warmup + last_restart &&
           rms > config.divergence_factor * best_rms
            if restarts >= config.max_restarts
                reason = :diverged
                options.verbose && @printf("INR %3d  RMS %.5g exceeds %.3g x the best %.5g, stopping\n",
                                           epoch, rms, config.divergence_factor, best_rms)
                break
            end
            restarts += 1
            last_restart = epoch
            rate /= 2
            ps        = deepcopy(best_ps)
            opt_state = Optimisers.setup(Optimisers.Adam(_inr_rate(config, rate, epoch, options.max_iter)), ps)
            options.verbose && @printf("INR %3d  RMS %.5g exceeds %.3g x the best %.5g, restarting at rate %.3g\n",
                                       epoch, rms, config.divergence_factor, best_rms, rate)
            continue
        end
        epoch == options.max_iter && break
        grad_m  = inv3d_gradient(problem, state)                 # analytic dL/dm = Jᵀr, G never traced
        grad_ps = back(Float32.(grad_m))[1]                      # dL/dtheta via the network's pullback only
        Optimisers.adjust!(opt_state, _inr_rate(config, rate, epoch, options.max_iter))
        opt_state, ps = Optimisers.update(opt_state, ps, grad_ps)
    end

    m_final     = best_m
    final_state = _inv3d_evaluate(problem, m_final, 0.0)
    rho         = copy(model.A)
    rho[cells]  = m_final
    pred        = _inv3d_J(problem, m_final) .* sigma .+ dfixed
    Inv3DResult(rho, pred, final_state.phid, sqrt(final_state.phid / nd), 0.0, ones(length(cells)), history,
                cells, reason == :target_rms, reason, config)
end

#---------- file workflow ----------

# shared options and the algorithm config from an InvCtrl.INR file
function _inr_from_ctrl(c::InvCtrlINR3D)
    options = Inv3DOptions(max_iter = c.max_iter, target_rms = c.target_rms, bounds = c.bounds)
    config = INR3DConfig(width = c.width, layers = c.layers, frequencies = c.frequencies,
                         learning_rate = c.learning_rate, scaling = c.scaling, seed = c.seed,
                         divergence_factor = c.divergence_factor, max_restarts = c.max_restarts,
                         schedule = c.schedule, rate_floor = c.rate_floor, anneal = c.anneal)
    options, config
end
