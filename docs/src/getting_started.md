# Getting Started

This page takes you from a fresh install to a density model, using one of the synthetic
surveys shipped with the package. Everything here runs in a few minutes on a laptop.

## Install

```julia
pkg> add GravityGeophysics
```

GravityGeophysics.jl needs **Julia 1.13.0**. If you are working from a clone of the
repository instead:

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

## The shape of a job

A gravity inversion here needs six files, all of them small text files sitting in one folder
next to your data:

| File | What it is |
|---|---|
| `data.xyz` | your survey: position, elevation, gravity, error |
| `model.start` | the mesh, filled with a starting density contrast |
| `model.prior` | what the model is pulled towards when the data say nothing |
| `FwdCtrl` | forward-solver settings |
| `InvCtrl.GN` or `InvCtrl.INR` | which inversion to run and how |
| `C3.dat` | which cells are free to change |

You do not have to write these by hand. [`MakeMesh3D`](data/mesh.md) builds the mesh, start
model and `C3.dat` from your survey, and the control files are short enough to copy from an
example and edit.

## Run the shipped example

From a clone of the repository:

```bash
julia --project=. -t auto examples/demo1/00_run_inv_GN.jl
julia --project=. -t auto examples/demo1/00_run_inv_INR.jl
julia --project=. examples/demo1/02_plot_model_XY.jl
```

The first two invert the same synthetic survey, once with Gauss–Newton and once with the
neural inversion, each into its own timestamped run folder. The third plots the true model
against both results; running and plotting are separate steps, so you can inspect a run before
deciding what to plot. `-t auto` lets Julia use all your cores, which is where most of the
speed comes from.

## Doing it yourself, step by step

```julia
using GravityGeophysics

# 1. design a mesh around the survey; writes model.start and C3.dat next to the data
MakeMesh3D("data.xyz")

# 2. load what you need
model = ReadModel3D("model.start")
data  = ReadData3D("data.xyz")
fwd   = FwdCtrl3D()

# 3. predict the gravity of a model
pred = forward3d(model, data)

# 4. invert, Gauss-Newton
gn = Invert3D(model, data, fwd;
              options = Inv3DOptions(max_iter = 20, target_rms = 1.0,
                                     bounds = (-500.0, 500.0)))

# 5. or invert with the neural representation; max_iter counts training epochs,
#    so it wants thousands rather than tens
inr = INR3D(model, data, fwd;
            config  = INR3DConfig(width = 256, layers = 4),
            options = Inv3DOptions(max_iter = 3000, target_rms = 1.0,
                                   bounds = (-500.0, 500.0)))
```

## The production path

In practice you will use the six-file call, which keeps a record of everything it did:

```julia
run = Invert3D("model.start", "data.xyz", "FwdCtrl",
               "InvCtrl.GN", "C3.dat", "model.prior")

PlotModelXY(run.final, "data.xyz")
```

This writes a **timestamped** folder next to the data, so a second run never overwrites the
first:

```
run_20261026_142233/
├── model.den        the recovered density contrast
├── data.pred        what that model predicts
├── History.log      rms and other diagnostics per iteration, as a fixed-width text table
├── Summary.txt      how it finished and under what settings
└── inputs/          copies of the six files it was given
```

Plotting is not part of the run: see [Plots](plots.md) for `PlotModelXY`/`XZ`/`YZ`/`XYZ`,
`PlotData3D` and `plot_inv3d_convergence`, the way `examples/run_inv3D.jl` and the demo
`01_plot_data.jl`/`0N_plot_model_*.jl` scripts write their own `plots/` folder beside the run.

To run the neural inversion instead, pass an `InvCtrl.INR` where the `InvCtrl.GN` was.
Nothing else changes — `ReadInvCtrl` picks the algorithm from the `Algorithm` line inside the
file.

## Did it work?

Three quick checks, in order:

1. **`Summary.txt`** — did it say `Converged: true`, and what RMS did it reach? An RMS near 1
   means the model explains the data to about the size of the errors.
2. **The data-fit figure** — the residual map should look like noise. Structure left in it
   means the model has not explained something.
3. **The convergence curve** — RMS should fall and flatten. If it is still falling steeply at
   the last iteration, give it more iterations.

## Next steps

- [Survey and model files](data/files.md) — the formats, in detail
- [Building a mesh](data/mesh.md) — what `MakeMesh3D` decides for you, and how to override it
- [Control files](data/control_files.md) — every setting, with sensible defaults
- [Choosing between the inversions](inversion/choosing.md)
