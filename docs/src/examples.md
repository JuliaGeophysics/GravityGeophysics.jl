# Worked examples

The `examples/` folder holds two complete synthetic surveys and a set of short scripts for the
everyday jobs.

## A note on the data

Both surveys are **synthetic**. Every model, station position and gravity value in them is
computed, not measured. They are placed on the Finnish national grid (ETRS-TM35FIN,
`EPSG:3067`) around the GTK office in Otaniemi, Espoo, so that coordinates, ground elevations,
anomaly amplitudes and error bars read like a real ground survey in Finnish Precambrian
bedrock. No measured data is distributed with this package.

## Run them

Each demo folder splits running from plotting. `00_run_inv_GN.jl` and `00_run_inv_INR.jl` invert
the survey into a timestamped `run_gn_<timestamp>/` or `run_inr_<timestamp>/`; nothing is
plotted automatically:

```bash
julia --project=. -t auto examples/demo1/00_run_inv_GN.jl
julia --project=. -t auto examples/demo1/00_run_inv_INR.jl
```

Then plot whatever runs you have, on demand:

```bash
julia --project=. examples/demo1/01_plot_data.jl
julia --project=. examples/demo1/02_plot_model_XY.jl
julia --project=. examples/demo1/03_plot_model_XZ.jl
julia --project=. examples/demo1/04_plot_model_YZ.jl
```

Each plotting script finds the latest GN and INR run next to it and compares them against
`model.true`; run only one algorithm and it plots just that one. The run and plot folders are
gitignored, so you always reproduce them locally rather than reading someone else's output.

## `demo1` — blocky targets, dense grid

289 stations on a 100 m grid over 1.6 × 1.6 km of gently rolling coastal bedrock (ground
4–38 m). Two uniform ±200 kg/m³ blocks sit at 100–400 m depth: a mafic intrusion in granitic
gneiss, and a clay-filled bedrock depression. The anomaly runs −0.75 to +0.73 mGal with
2 % + 5 µGal errors.

This is the classical setting, and Gauss–Newton takes it: correlation 0.66 against the INR's
0.63, reaching RMS 0.99 in a single step where the INR needs 281 epochs.

Both overshoot the mass and put it too deep — 15.8 Mt at 447 m and 12.6 Mt at 403 m against a
true 8.8 Mt at 245 m. That is ordinary gravity depth ambiguity, not a failure of either
method, and it is worth seeing on a case where you know the answer.

## `demo2` — gradational target, sparse survey

A basin 280 m deep whose fill grades from −200 kg/m³ at the surface to background at the
floor, beside a 200 × 200 × 193 m block of +260 kg/m³ — a glacially scoured bedrock
depression filled with soft sediment next to a mafic body, about 2 km north-west of the
`demo1` site. 150 scattered stations, with no access inside a 130 m clearing between the
two bodies: 12 600 cells to constrain with 150 data.

Two things here are chosen to suit a coordinate network, and both are ordinary field
conditions rather than tricks: the truth is a **continuous function of position** rather than
a set of uniform blocks, and the survey is **sparse and irregular**.

Scored over the core against a true −3.05 Mt deficit at 52 m and +2.01 Mt excess at 147 m:

| | RMS | iterations | correlation | deficit | its centroid | excess | its centroid |
|---|---|---|---|---|---|---|---|
| Gauss–Newton | 0.95 | 3 | 0.61 | −6.34 Mt | 189 m | +2.52 Mt | 188 m |
| INR | 0.99 | 257 | **0.81** | **−2.89 Mt** | **54 m** | +1.52 Mt | **119 m** |

Gauss–Newton buys the misfit with twice the mass it needs, buried below 180 m, and speckles
the top layer where stations are sparse. The network returns a coherent bowl that shrinks with
depth like the true one, and a compact block in place.

Gauss–Newton is not strawmanned: across stronger smallness, heavier smoothing, depth weighting
and no weighting, its correlation stays between 0.41 and 0.65 and the basin centroid never
comes shallower than 149 m.

The data are forward-modelled on a 2× refinement of the inversion mesh and given 2 % + 3 µGal
noise, so the forward model sees structure no 40 m cell can hold — the true model *on the
inversion mesh* only fits the data at RMS 2.14. That is deliberate: a synthetic where the
truth fits perfectly is easier than any real survey.

## The short scripts

Each has a settings block at the top to edit, and sensible defaults so it runs as-is.

| Script | What it does |
|---|---|
| `make_mesh3D.jl` | design a mesh, start model and covariance around a survey |
| `run_fwd3D.jl` | predict gravity from a model at a survey's stations |
| `run_inv3D.jl` | run an inversion and plot it; the control file picks the algorithm |
| `plot_data3D.jl` | map a survey, or observed/predicted/residual together |
| `plot_model3D.jl` | open a model in a slice, section or 3D viewer |

```bash
julia --project=. examples/make_mesh3D.jl examples/demo1/data.xyz
julia --project=. -t auto examples/run_inv3D.jl examples/demo1 InvCtrl.INR
julia --project=. examples/plot_model3D.jl examples/demo1/model.true xz
```

## What ships, and what you make

Each demo folder holds **inputs only** — the files a six-file `Invert3D` call needs, plus
`model.true` so you can score the result yourself. Everything a run produces (run folders,
outputs, figures) is gitignored and regenerated when you run it.
