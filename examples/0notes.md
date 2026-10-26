# Examples

## Short scripts

Each has a settings block at the top and defaults that run as-is, so they double as templates
for your own work.

| Script | What it does |
|---|---|
| `make_mesh3D.jl` | design a mesh, start model and covariance around a survey |
| `run_fwd3D.jl` | predict gravity from a model at a survey's stations |
| `run_inv3D.jl` | run an inversion and plot it; the control file picks the algorithm |
| `plot_data3D.jl` | map a survey, or observed/predicted/residual together |
| `plot_model3D.jl` | open a model in a slice, section or 3D viewer (`xy`, `xz`, `yz`, `xyz`) |

```bash
julia --project=. examples/make_mesh3D.jl examples/demo1/data.xyz
julia --project=. -t auto examples/run_fwd3D.jl examples/demo1/model.true examples/demo1/data.xyz
julia --project=. -t auto examples/run_inv3D.jl examples/demo1 InvCtrl.INR
julia --project=. examples/plot_model3D.jl examples/demo1/model.true xz
```

## The two demo surveys

Two worked synthetics. Each folder is **inputs only** — the files a six-file `Invert3D` call
needs, plus the true model to score against. Running and plotting are separate scripts, so
nothing plots until you ask it to. Everything a run or plot produces (`run_gn_<timestamp>/`,
`run_inr_<timestamp>/`, their `inputs/`, `plots/` and output files) is gitignored and
reproduced locally:

```bash
julia --project=. -t auto examples/demo1/00_run_inv_GN.jl
julia --project=. -t auto examples/demo1/00_run_inv_INR.jl
julia --project=. examples/demo1/01_plot_data.jl
julia --project=. examples/demo1/02_plot_model_XY.jl
julia --project=. examples/demo1/03_plot_model_XZ.jl
julia --project=. examples/demo1/04_plot_model_YZ.jl
```

Every run is stamped with its start time, so repeated runs accumulate side by side instead of
overwriting one another. Each plotting script finds the most recent `run_gn_*`/`run_inr_*`
folder next to it; run only one algorithm and the plots compare `model.true` against just that
one.

## Files in a demo folder

| File | Role |
|---|---|
| `data.xyz` | the survey: x, y, elevation, gz (mGal), error, lat, lon |
| `model.start` | `.den` start model, the mesh the inversion runs on |
| `model.prior` | reference model for the smallness term |
| `model.true` | the synthetic truth, for scoring and for outlining bodies in the plots |
| `FwdCtrl` | forward control: sensitivity storage, Frechet output |
| `InvCtrl.GN` | Gauss–Newton control |
| `InvCtrl.INR` | INR control |
| `C3.dat` | model covariance, written by `MakeMesh3D` |
| `00_run_inv_GN.jl` / `00_run_inv_INR.jl` | invert the survey, no plotting |
| `01_plot_data.jl` | map the survey and the latest runs' data fit and convergence |
| `02_plot_model_XY.jl` / `03_plot_model_XZ.jl` / `04_plot_model_YZ.jl` | compare True/GN/INR on one slice or section |

Both algorithms read the same `C3.dat`: the covariance describes the mesh, not the algorithm.
Whether the lateral padding is free is the algorithm's own choice and lives in its control
file, as `Invert padding`. GN sets it `yes` — its sensitivity weighting prices the volume of a
padding cell, so those cells cost what they are worth. The INR sets it `no`, the default for
`InvCtrl.INR`: nothing in the INR pays for that volume, so with the padding free the network
fills those oversized far cells to buy misfit and the run stalls. `Summary.txt` records which
was used. The core is found from the mesh itself (`lateral_core_ranges`), the same ranges the
plotting and scoring already use.

## The demos

Both are **synthetic** — every model, station position and gravity value here is computed, not
measured. They sit on the Finnish national grid (ETRS-TM35FIN, `EPSG:3067`) around the GTK
office in Otaniemi, Espoo, so coordinates, elevations, amplitudes and errors read like a real
ground survey in Finnish Precambrian bedrock. No measured data ships with this package.

- **`demo1`** — 289 stations on a 100 m grid over 1.6 × 1.6 km of rolling coastal
  bedrock. Two uniform ±200 kg/m³ blocks at 100–400 m depth (a mafic intrusion in granitic
  gneiss, and a clay-filled bedrock depression), anomaly −0.75 to +0.73 mGal, errors
  2 % + 5 µGal. Blocky truth and dense data, the case that suits GN.
- **`demo2`** — a graded basin (−200 kg/m³ at surface, grading to background at the
  280 m floor) beside a +260 kg/m³ block, 150 scattered stations with a 130 m no-access
  clearing, ~2 km north-west of the `demo1` site. A truth that is a continuous function
  of position, sampled sparsely: the case that suits the INR, which wins on correlation and
  recovered mass.
