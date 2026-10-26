# 3D control files
# Author: @pankajkmishra
# FwdCtrl, InvCtrl (GN or INR) and the model covariance of the 3D gravity workflows
# Control files are "Key : value" lines; unknown, duplicate and missing required keys are errors

using Printf

#---------- key-value parsing ----------

_ctrl_key(s) = lowercase(join(split(strip(s)), " "))

function _ctrl_parse_bool(s)
    v = lowercase(strip(s))
    v in ("yes", "true", "1") && return true
    v in ("no", "false", "0") && return false
    throw(ArgumentError("expected yes or no, got '$s'"))
end

function _ctrl_parse_pair(s)
    t = split(strip(s))
    length(t) == 2 || throw(ArgumentError("expected two numbers, got '$s'"))
    (parse(Float64, t[1]), parse(Float64, t[2]))
end

# spec: display key => (field, parser, required)
function _read_ctrl(path::AbstractString, spec, kind::AbstractString)
    isfile(path) || error("$kind file not found: $path")
    lookup = Dict(_ctrl_key(k) => (k, v...) for (k, v) in spec)
    values = Dict{Symbol, Any}()
    for (n, raw) in enumerate(eachline(path))
        line = strip(first(split(raw, '#')))
        isempty(line) && continue
        occursin(':', line) || error("$path:$n: expected 'Key : value', got '$line'")
        key, value = strip.(split(line, ':'; limit = 2))
        entry = get(lookup, _ctrl_key(key), nothing)
        entry === nothing && error("$path:$n: unknown $kind key '$key'")
        _, field, parser, _ = entry
        haskey(values, field) && error("$path:$n: duplicate $kind key '$key'")
        values[field] = try
            parser(value)
        catch err
            error("$path:$n: bad value for '$key': $(sprint(showerror, err))")
        end
    end
    for (k, (field, _, required)) in spec
        required && !haskey(values, field) && error("$path: missing required $kind key '$k'")
    end
    values
end

function _write_ctrl(io, rows)
    width = maximum(length ∘ first, rows)
    foreach(((k, v),) -> println(io, rpad(k, width), " : ", v), rows)
end

_ctrl_float(s) = parse(Float64, s)
_ctrl_int(s) = parse(Int, s)
_ctrl_yesno(b::Bool) = b ? "yes" : "no"
_ctrl_g(x) = @sprintf("%.6g", x)

_ctrl_algorithm(s) = (a = Symbol(lowercase(strip(s))); a in (:gn, :inr) ? a :
                      throw(ArgumentError("algorithm must be GN or INR")))
_ctrl_lambda(s) = lowercase(strip(s)) == "auto" ? nothing : parse(Float64, s)
_ctrl_weighting(s) = (w = Symbol(lowercase(strip(s))); w in (:sensitivity, :depth, :none) ? w :
                      throw(ArgumentError("sensitivity weighting must be sensitivity, depth or none")))
_ctrl_depth0(s) = lowercase(strip(s)) == "auto" ? nothing : parse(Float64, s)
_ctrl_scaling(s) = (c = Symbol(replace(lowercase(strip(s)), ' ' => '_', '-' => '_')); c in (:per_axis, :isotropic) ? c :
                    throw(ArgumentError("INR coordinate scaling must be per-axis or isotropic")))
_ctrl_schedule(s) = (c = Symbol(lowercase(strip(s))); c in (:cosine, :none) ? c :
                     throw(ArgumentError("INR learning rate schedule must be cosine or none")))

#---------- FwdCtrl ----------

"""
    FwdCtrl3D(; write_frechet=false) -> FwdCtrl3D

Input of `forward3d()` and of every inversion: whether the sensitivity matrix is kept after the run.
"""
Base.@kwdef struct FwdCtrl3D
    write_frechet::Bool = false          # keep G as sensitivity.jld2 (Float32) after the run
end

const _FWD_CTRL_SPEC = (
    "Write Frechet derivative"   => (:write_frechet, _ctrl_parse_bool, false),
)

_validate_ctrl(c::FwdCtrl3D) = c

"""
    ReadFwdCtrl3D(path) -> FwdCtrl3D

Read a 3D forward control file.
"""
ReadFwdCtrl3D(path::AbstractString) = _validate_ctrl(FwdCtrl3D(; _read_ctrl(path, _FWD_CTRL_SPEC, "FwdCtrl")...))

"""
    WriteFwdCtrl3D(path, ctrl::FwdCtrl3D) -> path

Write a 3D forward control file.
"""
function WriteFwdCtrl3D(path::AbstractString, c::FwdCtrl3D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    rows = [("Write Frechet derivative", _ctrl_yesno(c.write_frechet))]
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- InvCtrl.GN ----------

"""
    InvCtrl3D(; algorithm, lambda, target_rms, max_iter, lambda_cooling=2.0,
                bounds=(-1000.0, 1000.0), max_linesearch=12, smallness=0.01,
                smooth_x=1.0, smooth_y=1.0, smooth_z=1.0, weighting=:sensitivity,
                depth_exponent=2.0, depth_z0=nothing, error_floor=0.0, error_percent=0.0,
                gn_damping=1e-3, cg_iterations=50, cg_tolerance=1e-3, invert_padding=true)

Input of a Gauss–Newton run, the in-memory form of an `InvCtrl.GN` file.
"""
Base.@kwdef struct InvCtrl3D
    algorithm::Symbol                                 # :gn
    lambda::Union{Nothing, Float64}                   # starting β; nothing ("auto") picks it from the problem
    target_rms::Float64                               # stop once the rms misfit falls below this
    max_iter::Int                                     # maximum accepted iterations
    lambda_cooling::Float64 = 2.0                     # β is divided by this per accepted step; 1 holds it fixed
    bounds::Tuple{Float64, Float64} = (-1000.0, 1000.0) # density contrast limits (kg/m³), enforced by projection
    max_linesearch::Int     = 12                      # backtracking halvings allowed per direction
    smallness::Float64      = 0.01                    # pull towards the prior; 1/L² sets a length of L core cells
    smooth_x::Float64       = 1.0                     # north-south smoothing weight, relative to a core cell
    smooth_y::Float64       = 1.0                     # east-west smoothing weight
    smooth_z::Float64       = 1.0                     # vertical smoothing weight; raise it against depth streaks
    weighting::Symbol       = :sensitivity            # :sensitivity, :depth or :none
    depth_exponent::Float64 = 2.0                     # β of the depth weighting (z + z₀)^(-β/2)
    depth_z0::Union{Nothing, Float64} = nothing       # z₀ of that weighting; nothing = median sensor height + half a cell
    error_floor::Float64    = 0.0                     # smallest error believed (mGal)
    error_percent::Float64  = 0.0                     # error as a percentage of |datum|; the largest of the three wins
    gn_damping::Float64     = 1e-3                    # Levenberg-Marquardt damping relative to the diagonal
    cg_iterations::Int      = 50                      # maximum inner conjugate-gradient iterations per step
    cg_tolerance::Float64   = 1e-3                    # inner conjugate-gradient tolerance
    invert_padding::Bool    = true                    # false fixes the lateral padding on top of the covariance mask
end

const _INV_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm, true),
    "Initial damping factor lambda"      => (:lambda, _ctrl_lambda, true),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Lambda cooling factor"              => (:lambda_cooling, _ctrl_float, false),
    "Density bounds (kg/m3)"             => (:bounds, _ctrl_parse_pair, false),
    "Max line search steps"              => (:max_linesearch, _ctrl_int, false),
    "Smallness weight"                   => (:smallness, _ctrl_float, false),
    "Smoothing weight x"                 => (:smooth_x, _ctrl_float, false),
    "Smoothing weight y"                 => (:smooth_y, _ctrl_float, false),
    "Smoothing weight z"                 => (:smooth_z, _ctrl_float, false),
    "Sensitivity weighting"              => (:weighting, _ctrl_weighting, false),
    "Depth weighting exponent"           => (:depth_exponent, _ctrl_float, false),
    "Depth weighting z0 (m)"             => (:depth_z0, _ctrl_depth0, false),
    "Error floor (mGal)"                 => (:error_floor, _ctrl_float, false),
    "Error percent"                      => (:error_percent, _ctrl_float, false),
    "GN damping"                         => (:gn_damping, _ctrl_float, false),
    "Maximum CG iterations"              => (:cg_iterations, _ctrl_int, false),
    "CG tolerance"                       => (:cg_tolerance, _ctrl_float, false),
    "Invert padding"                     => (:invert_padding, _ctrl_parse_bool, false),
)

function _validate_ctrl(c::InvCtrl3D)
    c.algorithm == :gn || throw(ArgumentError("InvCtrl.GN needs 'Algorithm : GN'; INR runs read InvCtrl.INR"))
    c.lambda === nothing || (isfinite(c.lambda) && c.lambda >= 0) || throw(ArgumentError("lambda must be auto or nonnegative"))
    c.lambda_cooling >= 1 || throw(ArgumentError("lambda cooling factor must be at least 1"))
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.target_rms >= 0 || throw(ArgumentError("target rms must be nonnegative"))
    c.bounds[1] < c.bounds[2] || throw(ArgumentError("density bounds must be increasing"))
    c.max_linesearch > 0 || throw(ArgumentError("max line search steps must be positive"))
    all(w -> isfinite(w) && w >= 0, (c.smallness, c.smooth_x, c.smooth_y, c.smooth_z)) ||
        throw(ArgumentError("regularization weights must be finite and nonnegative"))
    c.smallness + c.smooth_x + c.smooth_y + c.smooth_z > 0 || throw(ArgumentError("at least one regularization weight must be positive"))
    c.depth_exponent >= 0 || throw(ArgumentError("depth weighting exponent must be nonnegative"))
    c.depth_z0 === nothing || c.depth_z0 >= 0 || throw(ArgumentError("depth weighting z0 must be auto or nonnegative"))
    c.error_floor >= 0 && c.error_percent >= 0 || throw(ArgumentError("error floor and percent must be nonnegative"))
    c.gn_damping >= 0 || throw(ArgumentError("GN damping must be nonnegative"))
    c.cg_iterations > 0 && c.cg_tolerance > 0 || throw(ArgumentError("CG iterations and tolerance must be positive"))
    c
end

"""
    ReadInvCtrl3D(path) -> InvCtrl3D

Read a GN inversion control file (`InvCtrl.GN`).
"""
ReadInvCtrl3D(path::AbstractString) = _validate_ctrl(InvCtrl3D(; _read_ctrl(path, _INV_CTRL_SPEC, "InvCtrl.GN")...))

"""
    WriteInvCtrl3D(path, ctrl::InvCtrl3D) -> path

Write a GN inversion control file with every key.
"""
function WriteInvCtrl3D(path::AbstractString, c::InvCtrl3D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g = _ctrl_g
    rows = [
        ("Algorithm", "GN"),
        ("Initial damping factor lambda", c.lambda === nothing ? "auto" : g(c.lambda)),
        ("Exit search when rms is less than", g(c.target_rms)),
        ("Maximum number of iterations", string(c.max_iter)),
        ("Lambda cooling factor", g(c.lambda_cooling)),
        ("Density bounds (kg/m3)", "$(g(c.bounds[1])) $(g(c.bounds[2]))"),
        ("Max line search steps", string(c.max_linesearch)),
        ("Smallness weight", g(c.smallness)),
        ("Smoothing weight x", g(c.smooth_x)),
        ("Smoothing weight y", g(c.smooth_y)),
        ("Smoothing weight z", g(c.smooth_z)),
        ("Sensitivity weighting", string(c.weighting)),
        ("Depth weighting exponent", g(c.depth_exponent)),
        ("Depth weighting z0 (m)", c.depth_z0 === nothing ? "auto" : g(c.depth_z0)),
        ("Error floor (mGal)", g(c.error_floor)),
        ("Error percent", g(c.error_percent)),
        ("GN damping", g(c.gn_damping)),
        ("Maximum CG iterations", string(c.cg_iterations)),
        ("CG tolerance", g(c.cg_tolerance)),
        ("Invert padding", _ctrl_yesno(c.invert_padding)),
    ]
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

#---------- InvCtrl.INR ----------

"""
    InvCtrlINR3D(; algorithm, target_rms, max_iter, bounds=(-1000.0, 1000.0),
                   error_floor=0.0, error_percent=0.0, width=256, layers=4, frequencies=2,
                   learning_rate=1e-3, scaling=:per_axis, seed=42, divergence_factor=3.0,
                   max_restarts=3, schedule=:cosine, rate_floor=0.05, anneal=0.5,
                   invert_padding=false)

Input of an INR run, the in-memory form of an `InvCtrl.INR` file.
"""
Base.@kwdef struct InvCtrlINR3D
    algorithm::Symbol                                 # :inr
    target_rms::Float64                               # stop as soon as the rms reaches this
    max_iter::Int                                     # maximum training epochs, not iterations
    bounds::Tuple{Float64, Float64} = (-1000.0, 1000.0) # density contrast limits (kg/m³), the network output range
    error_floor::Float64    = 0.0                     # smallest error believed (mGal)
    error_percent::Float64  = 0.0                     # error as a percentage of |datum|; the largest of the three wins
    width::Int              = 256                     # neurons per hidden layer; capacity, and so the regularization
    layers::Int             = 4                       # number of hidden layers
    frequencies::Int        = 2                       # Fourier feature octaves; more = sharper structure available
    learning_rate::Float64  = 1e-3                    # peak Adam step size
    scaling::Symbol         = :per_axis               # coordinate normalization, :per_axis or :isotropic
    seed::Int               = 42                      # network initialization; change it to test robustness
    divergence_factor::Float64 = 3.0                  # restart when the rms exceeds this multiple of the best; 0 disables
    max_restarts::Int       = 3                       # how often a diverging run may fall back and halve the rate
    schedule::Symbol        = :cosine                 # learning rate decay, :cosine or :none
    rate_floor::Float64     = 0.05                    # cosine decays the rate to this fraction of the peak
    anneal::Float64         = 0.5                     # fraction of the budget over which Fourier bands are released
    invert_padding::Bool    = false                   # the INR cannot pay for padding volume, so it is held fixed
end

const _INR_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm, true),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Density bounds (kg/m3)"             => (:bounds, _ctrl_parse_pair, false),
    "Error floor (mGal)"                 => (:error_floor, _ctrl_float, false),
    "Error percent"                      => (:error_percent, _ctrl_float, false),
    "Hidden width"                       => (:width, _ctrl_int, false),
    "Hidden layers"                      => (:layers, _ctrl_int, false),
    "Fourier frequencies"                => (:frequencies, _ctrl_int, false),
    "Learning rate"                      => (:learning_rate, _ctrl_float, false),
    "Coordinate scaling"                 => (:scaling, _ctrl_scaling, false),
    "Seed"                               => (:seed, _ctrl_int, false),
    "Divergence factor"                  => (:divergence_factor, _ctrl_float, false),
    "Maximum restarts"                   => (:max_restarts, _ctrl_int, false),
    "Learning rate schedule"             => (:schedule, _ctrl_schedule, false),
    "Learning rate floor"                => (:rate_floor, _ctrl_float, false),
    "Frequency annealing"                => (:anneal, _ctrl_float, false),
    "Invert padding"                     => (:invert_padding, _ctrl_parse_bool, false),
)

function _validate_ctrl(c::InvCtrlINR3D)
    c.algorithm == :inr || throw(ArgumentError("InvCtrl.INR needs 'Algorithm : INR'; GN runs read InvCtrl.GN"))
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c.target_rms >= 0 || throw(ArgumentError("target rms must be nonnegative"))
    c.bounds[1] < c.bounds[2] || throw(ArgumentError("density bounds must be increasing"))
    c.error_floor >= 0 && c.error_percent >= 0 || throw(ArgumentError("error floor and percent must be nonnegative"))
    c.width > 0 && c.layers > 0 && c.frequencies >= 0 ||
        throw(ArgumentError("hidden width and layers must be positive and the frequencies nonnegative"))
    c.learning_rate > 0 || throw(ArgumentError("learning rate must be positive"))
    c.divergence_factor == 0 || c.divergence_factor > 1 ||
        throw(ArgumentError("divergence factor must be zero or greater than one"))
    c.max_restarts >= 0 || throw(ArgumentError("the number of restarts must be nonnegative"))
    c.schedule in (:cosine, :none) || throw(ArgumentError("the schedule must be cosine or none"))
    0 < c.rate_floor <= 1 || throw(ArgumentError("the learning rate floor must be in (0, 1]"))
    0 <= c.anneal <= 1 || throw(ArgumentError("the frequency annealing must be in [0, 1]"))
    c
end

"""
    ReadInvCtrlINR3D(path) -> InvCtrlINR3D

Read an INR inversion control file (`InvCtrl.INR`).
"""
ReadInvCtrlINR3D(path::AbstractString) = _validate_ctrl(InvCtrlINR3D(; _read_ctrl(path, _INR_CTRL_SPEC, "InvCtrl.INR")...))

"""
    WriteInvCtrlINR3D(path, ctrl::InvCtrlINR3D) -> path

Write an INR inversion control file with every key.
"""
function WriteInvCtrlINR3D(path::AbstractString, c::InvCtrlINR3D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g = _ctrl_g
    rows = [
        ("Algorithm", "INR"),
        ("Exit search when rms is less than", g(c.target_rms)),
        ("Maximum number of iterations", string(c.max_iter)),
        ("Density bounds (kg/m3)", "$(g(c.bounds[1])) $(g(c.bounds[2]))"),
        ("Error floor (mGal)", g(c.error_floor)),
        ("Error percent", g(c.error_percent)),
        ("Hidden width", string(c.width)),
        ("Hidden layers", string(c.layers)),
        ("Fourier frequencies", string(c.frequencies)),
        ("Learning rate", g(c.learning_rate)),
        ("Coordinate scaling", replace(string(c.scaling), '_' => '-')),
        ("Seed", string(c.seed)),
        ("Divergence factor", g(c.divergence_factor)),
        ("Maximum restarts", string(c.max_restarts)),
        ("Learning rate schedule", string(c.schedule)),
        ("Learning rate floor", g(c.rate_floor)),
        ("Frequency annealing", g(c.anneal)),
        ("Invert padding", _ctrl_yesno(c.invert_padding)),
    ]
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end

"""
    ReadInvCtrl(path) -> InvCtrl3D or InvCtrlINR3D

Read an inversion control file of either algorithm, chosen by its `Algorithm` line.
"""
function ReadInvCtrl(path::AbstractString)
    isfile(path) || error("InvCtrl file not found: $path")
    for raw in eachline(path)
        line = strip(first(split(raw, '#')))
        occursin(':', line) || continue
        key, value = strip.(split(line, ':'; limit = 2))
        _ctrl_key(key) == "algorithm" || continue
        a = _ctrl_algorithm(value)
        return a == :gn ? ReadInvCtrl3D(path) : ReadInvCtrlINR3D(path)
    end
    error("$path: missing required InvCtrl key 'Algorithm'")
end

#---------- model covariance ----------

# the covariance reader skips exactly 16 header lines - keep this 16 lines long
const _MESH_COV_HEADER = """
+-----------------------------------------------------------------------------+
| This file defines model covariance for a recursive autoregression scheme.   |
| The model space may be divided into distinct areas using integer masks.     |
| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |
| air, ocean and the rest of the model is turned off automatically. You can   |
| also define exceptions to override smoothing between any two model areas.   |
| To turn off smoothing set it to zero. This header is 16 lines long.         |
| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |
| 2. Smoothing in the X direction (NzEarth real values)                       |
| 3. Smoothing in the Y direction (NzEarth real values)                       |
| 4. Vertical smoothing (1 real value)                                        |
| 5. Number of times the smoothing should be applied (1 integer >= 0)         |
| 6. Number of exceptions (1 integer >= 0)                                    |
| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 2 & 3) |
| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|
+-----------------------------------------------------------------------------+
"""
"""
    Cov3D(; sx, sy, sz, n_smooth, exceptions=[], mask)

Model covariance; the inversion reads only `mask`, its smoothing comes from `InvCtrl`.
"""
Base.@kwdef struct Cov3D
    sx::Vector{Float64}                        # lateral smoothing in x, one value per layer, top to bottom
    sy::Vector{Float64}                        # lateral smoothing in y, one value per layer
    sz::Float64                                # vertical smoothing
    n_smooth::Int                              # number of smoothing passes
    exceptions::Vector{Tuple{Int, Int, Float64}} = Tuple{Int, Int, Float64}[]  # (class_a, class_b, smoothing) rules
    mask::Array{Int, 3}                        # 0 = air or fixed, 1 = free, other values = exception classes
end
"""
    Cov3D(mask; smoothing=0.3, n_smooth=2)

Uniform covariance around a mask.
"""
Cov3D(mask::AbstractArray{<:Integer, 3}; smoothing::Real = 0.3, n_smooth::Integer = 2) =
    Cov3D(sx = fill(Float64(smoothing), size(mask, 3)), sy = fill(Float64(smoothing), size(mask, 3)),
          sz = Float64(smoothing), n_smooth = n_smooth, mask = Array{Int, 3}(mask))

function _validate_ctrl(c::Cov3D)
    nz = size(c.mask, 3)
    length(c.sx) == nz && length(c.sy) == nz || throw(ArgumentError("the covariance needs one x and y smoothing per layer"))
    all(s -> isfinite(s) && s >= 0, c.sx) && all(s -> isfinite(s) && s >= 0, c.sy) && isfinite(c.sz) && c.sz >= 0 || throw(ArgumentError("covariance smoothing values must be finite and nonnegative"))
    c.n_smooth >= 0 || throw(ArgumentError("the number of smoothing passes must be nonnegative"))
    all(>=(0), c.mask) || throw(ArgumentError("covariance mask values must be nonnegative"))
    c
end

"""
    ReadCov3D(path) -> Cov3D

Read a covariance file: the 16-line header, `nx ny nz`, `nz` x and `nz` y smoothings,
the z smoothing, the number of passes, the exception rules, then for each layer range two
layer indices and `nx` rows of `ny` integers.
"""
function ReadCov3D(path::AbstractString)
    isfile(path) || error("covariance file not found: $path")
    lines = readlines(path)
    length(lines) > 16 || error("$path: covariance file shorter than its 16-line header")
    tokens = String[]
    foreach(l -> append!(tokens, split(l)), lines[17:end])
    k = 0
    next() = (k += 1; k <= length(tokens) || error("$path: covariance file ends early"); tokens[k])
    nx, ny, nz = parse(Int, next()), parse(Int, next()), parse(Int, next())
    sx         = [parse(Float64, next()) for _ in 1:nz]
    sy         = [parse(Float64, next()) for _ in 1:nz]
    sz         = parse(Float64, next())
    n_smooth   = parse(Int, next())
    exceptions = [(parse(Int, next()), parse(Int, next()), parse(Float64, next())) for _ in 1:parse(Int, next())]
    mask       = fill(-1, nx, ny, nz)
    while k < length(tokens)
        k1, k2 = parse(Int, next()), parse(Int, next())
        1 <= k1 <= k2 <= nz || error("$path: bad layer range $k1 $k2")
        block = [parse(Int, next()) for _ in 1:nx*ny]
        layer = permutedims(reshape(block, ny, nx))      # nx rows of ny values
        for kk in k1:k2
            mask[:, :, kk] = layer
        end
    end
    any(<(0), mask) && error("$path: covariance mask does not cover every layer")
    _validate_ctrl(Cov3D(; sx, sy, sz, n_smooth, exceptions, mask))
end
"""
    WriteCov3D(path, cov::Cov3D) -> path

Write a covariance file, one block per layer.
"""
function WriteCov3D(path::AbstractString, c::Cov3D)
    _validate_ctrl(c)
    nx, ny, nz = size(c.mask)
    mkpath(dirname(abspath(path)))
    line(io, v) = (foreach(x -> @printf(io, " %.3f ", x), v); println(io))
    open(path, "w") do io
        print(io, _MESH_COV_HEADER)
        println(io)
        @printf(io, " %d %12d %12d \n\n", nx, ny, nz)
        line(io, c.sx)
        line(io, c.sy)
        line(io, (c.sz,))
        println(io)
        @printf(io, " %d \n\n", c.n_smooth)
        @printf(io, " %d \n", length(c.exceptions))
        for (a, b, v) in c.exceptions
            @printf(io, " %d %d %.3f\n", a, b, v)
        end
        println(io)
        for k in 1:nz
            @printf(io, "\n %d %12d \n", k, k)
            for ix in 1:nx
                for iy in 1:ny
                    @printf(io, " %d ", c.mask[ix, iy, k])
                end
                println(io)
            end
        end
    end
    String(path)
end
