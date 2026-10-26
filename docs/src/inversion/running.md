# Running an inversion

Inversion goes the other way from forward modelling: you have the gravity readings and want
the density model that produced them. This page is about the mechanics — which files, which
call, what comes out. Which *method* to use is [a separate question](choosing.md).

## The six-file call

This is the one to use for real work:

```julia
run = Invert3D("model.start", "data.xyz", "FwdCtrl",
               "InvCtrl.GN", "C3.dat", "model.prior")
```

The arguments are always in that order: start model, data, forward control, inversion control,
covariance, prior. To run the neural inversion instead, pass an `InvCtrl.INR` in place of the
`InvCtrl.GN`. Nothing else changes — the `Algorithm` line inside the file decides, and
everything downstream is identical.

## What you get back

Every run writes its own **timestamped** folder next to the data, so running again never
destroys what you had:

```
run_20261026_142233/
├── model.den       the recovered density contrast
├── data.pred       what that model predicts at your stations
├── History.log     rms, objective and step size per iteration, as a fixed-width text table
├── Summary.txt     termination reason, final rms, settings used
└── inputs/         copies of the six input files
```

Pass `run_dir = "trial_a"` to label a run and it becomes `trial_a_20261026_142233/`.

In Julia the call returns a named tuple:

```julia
run.run_dir        # where everything was written
run.final          # the recovered model
run.predicted      # predicted data
run.rms            # final rms misfit
run.converged      # did it reach the target?
run.reason         # :target_rms, :max_iter, :diverged
run.history        # per-iteration diagnostics
```

## Reading `Summary.txt`

```
Algorithm: GN
Termination: target_rms
Converged: true
RMS: 0.910920
Final lambda: 0.0954525
Model weighting: sensitivity
Data count: 289
Accepted iterations: 3
Active cells: 47876
Padding inverted: true
Air cells: 2124
Fixed earth cells: 0
Density contrast range: -32.95 to 23.79 kg/m3
```

What to look at:

- **`Termination`** — `target_rms` is what you want. `max_iter` means it ran out of iterations
  while still improving; give it more. `diverged` means something is wrong with the setup.
- **`RMS`** — roughly 1 means the model explains the data to within the error bars. Much below
  1 and you are fitting noise; much above and the model cannot explain the data.
- **`Active cells`** — how many cells were actually free. If this is far smaller than you
  expect, check your `C3.dat` mask and the `Invert padding` setting.
- **`Density contrast range`** — if this is pinned against your `Density bounds`, the data are
  pushing harder than the bounds allow; widen them or regularize more.

## From memory instead of files

When you are experimenting it is often easier to stay in Julia:

```julia
model = ReadModel3D("model.start")
data  = ReadData3D("data.xyz")
fwd   = FwdCtrl3D()
opts  = Inv3DOptions(max_iter = 20, target_rms = 1.0, bounds = (-500.0, 500.0))

gn  = Invert3D(model, data, fwd; options = opts)
inr = INR3D(model, data, fwd; config = INR3DConfig(width = 256, layers = 4),
            options = Inv3DOptions(max_iter = 3000, target_rms = 1.0,
                                   bounds = (-500.0, 500.0)))
```

These skip the run folder, so nothing is written unless you write it.

## Making it faster

- Start Julia with `-t auto`. The sensitivity calculation threads over the data and it is the
  bulk of the time.
- Coarsen the mesh. Gravity cannot resolve much below your station spacing anyway, so fine
  cells often cost time without buying resolution.
- Coarsen the mesh if memory is the constraint; the sensitivity matrix is built in memory.

## When it will not converge

| Symptom | Usual cause |
|---|---|
| RMS flattens well above 1 | errors too optimistic, or regularization too strong. Check `Error floor`. |
| RMS plunges far below 1 | errors too pessimistic — you are fitting noise. |
| model pinned at the bounds | bounds too tight, or too little regularization |
| speckled, isolated cells | not enough smoothing, or no sensitivity/depth weighting |
| everything just under the stations | no depth or sensitivity weighting (GN) |
| INR stalls and will not improve | padding left free — check `Invert padding : no` |
