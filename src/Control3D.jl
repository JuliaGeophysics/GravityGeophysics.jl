# 3D control files
# Author: @pankajkmishra
# FwdCtrl and the inversion control (LSQR) of the 3D workflows
# Control files are "Key : value" lines; unknown, duplicate and missing required keys are errors

using Printf

#---------- key-value parsing ----------

_ctrl_key(s) = lowercase(join(split(strip(s)), " "))

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
_ctrl_solver(s) = (v = Symbol(lowercase(strip(s))); v in (:dense, :matrixfree) ? v :
                   throw(ArgumentError("forward solver must be dense or matrixfree")))
_ctrl_algorithm(s) = (a = Symbol(lowercase(strip(s))); a == :lsqr ? a : throw(ArgumentError("algorithm must be LSQR")))
_ctrl_weighting(s) = (w = Symbol(lowercase(strip(s))); w in (:none, :depth, :sensitivity) ? w :
                      throw(ArgumentError("model weighting must be none, depth or sensitivity")))
_ctrl_g(x) = @sprintf("%.6g", x)
_ctrl_pair_text(p) = "$(_ctrl_g(p[1])) $(_ctrl_g(p[2]))"

#---------- FwdCtrl ----------

"""
    FwdCtrl3D

Forward control, read from `FwdCtrl`.
- `solver`: `:dense` stores the sensitivity matrix (fast repeated forwards); `:matrixfree`
  recomputes the kernels in every product (large airborne surveys); both give the same
  values within round-off
"""
Base.@kwdef struct FwdCtrl3D
    solver::Symbol = :dense
end

const _FWD_CTRL_SPEC = (
    "Forward solver" => (:solver, _ctrl_solver, false),
)

_validate_ctrl(c::FwdCtrl3D) = (c.solver in (:dense, :matrixfree) ||
                                throw(ArgumentError("forward solver must be dense or matrixfree")); c)

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
    open(io -> _write_ctrl(io, [("Forward solver", string(c.solver))]), path, "w")
    String(path)
end

#---------- inversion control ----------

"""
    InvCtrl3D

Deterministic inversion control, read from `InvCtrl.LSQR`. Each inversion method has its own
control file with only its own keys, named by its `Algorithm` line; a new method (or norm) adds its
own control rather than extending this one. The lambda, rms and
iteration keys are those of the MT inversion control.
- `algorithm`: `:lsqr` (regularized least squares, each step solved by LSQR)
- `lambda`: initial regularization weight, as a multiple of the ratio of the data and model
  Hessians at the start (density and data units differ from survey to survey)
- `target_rms`, `max_iter`: stopping controls
- `lambda_cooling`, `cooling_interval`: λ is divided by the factor every `cooling_interval`
  iterations until the rms reaches the target
- `target_tolerance`: a step that lands below `(1 − tolerance) × target_rms` is redone with a
  larger λ, so the run stops at the target instead of fitting the noise
- `max_step`, `max_linesearch`: cap on the largest density change per iteration (kg/m³),
  backtracking halvings per direction
- `smallness`, `smooth_x`, `smooth_y`, `smooth_z`: weights of the regularization terms
- `weighting`, `depth_exponent`: cell weights `:depth` (Li & Oldenburg 1998, depth^(−exponent); the
  default, with exponent 2 for gravity, which decays as distance⁻²), `:sensitivity` or `:none`
- `bounds`: density bounds (kg/m³), kept by projection
- `lsqr_max_iter`, `lsqr_tolerance`: LSQR iterations and tolerance per solve
"""
Base.@kwdef struct InvCtrl3D
    algorithm::Symbol
    lambda::Float64
    target_rms::Float64
    max_iter::Int
    lambda_cooling::Float64 = 2.0
    cooling_interval::Int = 1
    target_tolerance::Float64 = 0.1
    max_step::Float64 = Inf
    max_linesearch::Int = 20
    smallness::Float64 = 0.01
    smooth_x::Float64 = 1.0
    smooth_y::Float64 = 1.0
    smooth_z::Float64 = 1.0
    weighting::Symbol = :depth
    depth_exponent::Float64 = 2.0
    bounds::Tuple{Float64, Float64} = (-1000.0, 1000.0)
    lsqr_max_iter::Int = 200
    lsqr_tolerance::Float64 = 1e-4
end

const _INV_CTRL_SPEC = (
    "Algorithm"                          => (:algorithm, _ctrl_algorithm, true),
    "Initial damping factor lambda"      => (:lambda, _ctrl_float, true),
    "Exit search when rms is less than"  => (:target_rms, _ctrl_float, true),
    "Maximum number of iterations"       => (:max_iter, _ctrl_int, true),
    "Lambda cooling factor"              => (:lambda_cooling, _ctrl_float, false),
    "Lambda cooling interval"            => (:cooling_interval, _ctrl_int, false),
    "Target rms tolerance"               => (:target_tolerance, _ctrl_float, false),
    "Max density step (kg/m3)"       => (:max_step, _ctrl_float, false),
    "Max line search steps"              => (:max_linesearch, _ctrl_int, false),
    "Smallness weight"                   => (:smallness, _ctrl_float, false),
    "Smoothing weight x"                 => (:smooth_x, _ctrl_float, false),
    "Smoothing weight y"                 => (:smooth_y, _ctrl_float, false),
    "Smoothing weight z"                 => (:smooth_z, _ctrl_float, false),
    "Model weighting"                    => (:weighting, _ctrl_weighting, false),
    "Depth exponent"                     => (:depth_exponent, _ctrl_float, false),
    "Density bounds (kg/m3)"         => (:bounds, _ctrl_parse_pair, false),
    "Maximum LSQR iterations"            => (:lsqr_max_iter, _ctrl_int, false),
    "LSQR tolerance"                     => (:lsqr_tolerance, _ctrl_float, false),
)

function _validate_ctrl(c::InvCtrl3D)
    c.algorithm == :lsqr || throw(ArgumentError("algorithm must be LSQR"))
    isfinite(c.lambda) && c.lambda > 0 || throw(ArgumentError("lambda must be finite and positive"))
    c.max_iter >= 0 || throw(ArgumentError("maximum number of iterations must be nonnegative"))
    c
end

"""
    ReadInvCtrl3D(path) -> InvCtrl3D

Read an LSQR inversion control file.
"""
ReadInvCtrl3D(path::AbstractString) = _validate_ctrl(InvCtrl3D(; _read_ctrl(path, _INV_CTRL_SPEC, "inversion control")...))

"""
    WriteInvCtrl3D(path, ctrl::InvCtrl3D) -> path

Write an LSQR inversion control file.
"""
function WriteInvCtrl3D(path::AbstractString, c::InvCtrl3D)
    _validate_ctrl(c)
    mkpath(dirname(abspath(path)))
    g = _ctrl_g
    rows = [
        ("Algorithm", uppercase(string(c.algorithm))),
        ("Initial damping factor lambda", g(c.lambda)),
        ("Exit search when rms is less than", g(c.target_rms)),
        ("Maximum number of iterations", string(c.max_iter)),
        ("Lambda cooling factor", g(c.lambda_cooling)),
        ("Lambda cooling interval", string(c.cooling_interval)),
        ("Target rms tolerance", g(c.target_tolerance)),
        ("Max density step (kg/m3)", g(c.max_step)),
        ("Max line search steps", string(c.max_linesearch)),
        ("Smallness weight", g(c.smallness)),
        ("Smoothing weight x", g(c.smooth_x)),
        ("Smoothing weight y", g(c.smooth_y)),
        ("Smoothing weight z", g(c.smooth_z)),
        ("Model weighting", string(c.weighting)),
        ("Depth exponent", g(c.depth_exponent)),
        ("Density bounds (kg/m3)", _ctrl_pair_text(c.bounds)),
        ("Maximum LSQR iterations", string(c.lsqr_max_iter)),
        ("LSQR tolerance", g(c.lsqr_tolerance)),
    ]
    open(io -> _write_ctrl(io, rows), path, "w")
    String(path)
end
