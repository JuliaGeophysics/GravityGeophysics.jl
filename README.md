<h1 align="center">GravityGeophysics.jl</h1>

<p align="center"><em>A software repository for gravity geophysics research and applications.</em></p>

<p align="center">
	<a href="https://github.com/JuliaGeophysics/GravityGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

3D forward modelling and inversion of gravity data for density contrast. Everything is driven by small
plain-text files: a `.xyz` survey, block-layout `.den` models, `Key : value` control files, and a run
folder written next to the data.

Two inversions share one driver, one control-file workflow and one result type:

- **Gauss–Newton** — projected, damped, matrix-free CG, with smallness/smoothing regularization and
  sensitivity or depth weighting.
- **Implicit neural representation** — the density of each cell is the bounded output of a coordinate MLP
  ([Lux.jl](https://github.com/LuxDL/Lux.jl)), trained by Adam against the data misfit alone; the network's
  capacity and early stopping are the regularization.

## Requirements

- [Julia](https://julialang.org) 1.13.0 — the only supported version

## Installation

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

## Workflow

```julia
using GravityGeophysics

MakeMesh3D("data.xyz")                      # padded mesh + start model + C3.dat from the survey

model = ReadModel3D("model.start")          # block layout, kg/m³, air as 1e17
data  = ReadData3D("data.xyz")              # x/y/elevation/value/error, mGal
fwd   = FwdCtrl3D()

pred = forward3d(model, data, fwd)          # predicted gz, mGal

gn  = Invert3D(model, data, fwd;            # Gauss–Newton
               options = Inv3DOptions(max_iter = 20, target_rms = 1.0, bounds = (-500.0, 500.0)))

inr = INR3D(model, data, fwd;               # INR: max_iter is Adam epochs, so it wants thousands
            config = INR3DConfig(width = 256, layers = 4),
            options = Inv3DOptions(max_iter = 3000, target_rms = 1.0, bounds = (-500.0, 500.0)))
```

The production path is the six-file call, which writes a run folder next to the data:

```julia
run = Invert3D("model.start", "data.xyz", "FwdCtrl", "InvCtrl.GN", "C3.dat", "model.prior")
```

Every run gets its own timestamped folder, so nothing is overwritten: `run_dir` is a label, stamped into
`<run_dir>_YYYYmmdd_HHMMSS/`, and without one the folder is `run_YYYYmmdd_HHMMSS/` next to the data. It
holds `model.den`, `data.pred`, `History.log`, `Summary.txt` and `inputs/`.
Passing an `InvCtrl.INR` instead of an `InvCtrl.GN` runs the INR; `ReadInvCtrl` dispatches on the
`Algorithm` key, and everything else about the workflow is identical. Threads (`-t auto`) parallelise the
sensitivity over the data. There is no inversion-plotting API: `PlotModelXY`, `PlotModelXZ`, `PlotModelYZ`
and `PlotModelXYZ` open a model (or a `label => model` list, side by side on one colour scale) in a
slider-driven window with an Export Figure button, and `PlotData3D`/`plot_inv3d_convergence` cover the
data fit and the convergence; a script composes these the way it wants its figures, the way
`examples/run_inv3D.jl` and the demo scripts do.

## Examples

`examples/` holds two synthetic demos, inputs only: the nine files a six-file call needs, including
`model.true` to score the result against. Each folder splits running from plotting: `00_run_inv_GN.jl` and
`00_run_inv_INR.jl` invert the data into `run_gn_<timestamp>/` and `run_inr_<timestamp>/`, and
`01_plot_data.jl`, `02_plot_model_XY.jl`, `03_plot_model_XZ.jl` and `04_plot_model_YZ.jl` plot whatever runs
they find next to the latest GN/INR result, into `plots/`. Run folders, outputs and plots are gitignored and
reproduced locally; the scripts that built the synthetics live outside the repository.

```bash
julia --project=. -t auto examples/demo1/00_run_inv_GN.jl
julia --project=. -t auto examples/demo1/00_run_inv_INR.jl
julia --project=. examples/demo1/02_plot_model_XY.jl
```

Both are **synthetic**: the models, the station positions and every gravity value are computed, not
measured. They are placed on the Finnish national grid (ETRS-TM35FIN, `EPSG:3067`) around the GTK office in
Otaniemi, Espoo, so that coordinates, ground elevations, anomaly amplitudes and error bars read like a real
ground survey in Finnish Precambrian bedrock. No measured data is distributed with this package.

**`demo1` — two blocks under topography, the case that suits GN.** 289 stations on a 100 m grid over
1.6 × 1.6 km of gently rolling coastal bedrock (ground 4–38 m). Two uniform ±200 kg/m³ blocks sit at
100–400 m depth: a mafic intrusion against granitic gneiss, and a clay-filled bedrock depression. The
anomaly runs −0.75 to +0.73 mGal with 2 % + 5 µGal errors (7.9 µGal mean), which is what a careful ground
survey with a modern relative meter delivers. Blocky truth and dense data is the classical setting, and GN
takes it: correlation 0.66 against the INR's 0.63, reaching RMS 0.99 in a single step where the INR needs
281 epochs for RMS 1.00. Both overshoot the mass and put it too deep — GN returns a 15.8 Mt deficit at
447 m and the INR 12.6 Mt at 403 m, against the true 8.8 Mt at 245 m — which is ordinary gravity depth
ambiguity, not a failure of either method.

**`demo2` — a graded basin beside an intrusion, the case that suits the INR.** A smooth bowl 280 m deep
whose fill goes from −200 kg/m³ at the surface to the background at the floor, next to a 200 × 200 × 193 m
block of +260 kg/m³, on flat ground — a glacially scoured bedrock depression filled with soft sediment
beside a mafic body, about 2 km north-west of the `demo1` site. Two things are chosen to favour a
coordinate network, and both are ordinary field conditions rather than tricks:

- the truth is a **continuous function of position**, not a set of uniform blocks — what an INR represents
  naturally and what `demo1` does not test;
- the survey is **sparse and irregular**: 150 scattered ground stations, not a filled grid, with no access
  inside a 130 m clearing between the two bodies — 12 600 cells to constrain with 150 data.

The mesh is 30 × 30 × 14 with an 18 × 18 core of 40 m cells. The data are forward-modelled on a 2×
refinement of that mesh and given 2 % + 3 µGal noise, so the forward model sees structure no 40 m cell can
hold: the true model *on the inversion mesh* itself only fits them at RMS 2.14.

Both algorithms reach the target. Scored over the core against the true −3.05 Mt deficit at a 52 m centroid
and +2.01 Mt excess at 147 m:

| | RMS | iterations | correlation | deficit | its centroid | excess | its centroid |
|---|---|---|---|---|---|---|---|
| GN (`InvCtrl.GN`, sensitivity weighting) | 0.95 | 3 | 0.61 | −6.34 Mt | 189 m | +2.52 Mt | 188 m |
| INR (`InvCtrl.INR`, 128 × 3, 4 frequencies) | 0.99 | 257 | **0.81** | **−2.89 Mt** | **54 m** | +1.52 Mt | **119 m** |

GN buys the misfit with twice the mass it needs, buried below 180 m, and speckles the top layer with
isolated cells at the bounds where stations are sparse (its `plots/ModelXY.png`); the network returns
a coherent bowl that shrinks with depth like the true one and a compact block in place, slightly light on
both. GN is not strawmanned here: across stronger smallness, heavier smoothing, depth weighting and no
weighting its correlation stays between 0.41 and 0.65, and the basin centroid never comes shallower than
149 m against the true 52 m.

Both runs read the same `C3.dat`; the INR's control file holds the padding fixed with
`Invert padding : no`.

## Conventions and files

- **Notation**: Blakely (1995) throughout — `γ` the gravitational constant (`GAMMA` in the code, leaving `G`
  for the sensitivity matrix), `ρ` the density contrast, `Φ` the gravitational potential, `g` the field.
- **Frame**: x north, y east, depth positive down from the top of the mesh (the highest topography),
  whose elevation is the `# Top elevation (m)` header of the model file. A sensor's depth is
  `top_elevation - elevation`.
- **Data** (`.xyz`): `easting northing elevation value error`, then latitude/longitude (WGS84) columns
  derived from the CRS when absent, and `# Key : value` header lines (`CRS`, `Geographic CRS`, title).
  Geographic coordinates are projected to UTM on reading.
- **Model** (`.den`): the block layout, `nx ny nz 0 LINEAR`, density contrast in kg/m³, air above the
  topography as 1e17 on file and NaN in memory.
- **FwdCtrl**: `Write Frechet derivative`, which keeps the sensitivity as `sensitivity.jld2` in the run folder.
- **InvCtrl.GN**: lambda, rms and iteration keys, plus density bounds, smallness and smoothing
  weights, sensitivity or depth weighting, error floor and percentage, β cooling and the CG solve.
- **InvCtrl.INR**: target rms, epochs, bounds, error floor and percentage, plus the network's width, layers,
  Fourier frequencies, peak learning rate, coordinate scaling and seed, and the training schedule
  (`Learning rate schedule`, `Learning rate floor`, `Frequency annealing`). A file that leaves the three
  training keys out reads with their defaults, so control files written before 0.3 still load.
- **C3.dat**: the model covariance; the inversion reads only its mask (0 = air or fixed, otherwise free). It
  describes the mesh, so one file serves both algorithms — whether the lateral padding is inverted is the
  algorithm's choice and lives in its control file as `Invert padding`.

## Method

The forward solver sums the closed-form gz of uniform rectangular prisms (Plouff 1976; Blakely 1995;
Nagy et al. 2000), evaluated once per mesh node and differenced over the corners of every cell; a
finite-difference Poisson alternative (`forward3d_fd`) solves `∇²Φ = 4πγρ` with a node-centred
finite-volume operator, monopole Dirichlet boundaries and CG with Jacobi preconditioning. The model's grid
is padded with air on every face, so the sensors — which sit above the mesh, not in it — and their
differencing stencils lie inside the solved domain; the operator is symmetric on any tensor grid, so
stretched and laterally padded meshes and topography (air cells carry no mass) are taken as they are. It
agrees with the prism sum to about 1–2 % RMS on uniform meshes and 8 % on `demo2`'s stretched one.

The Gauss–Newton inversion minimises `½‖Wd(Gm − d)‖² + ½β‖R W (m − mref)‖²` under density bounds
(Li & Oldenburg 1996), with sensitivity (Li & Oldenburg 2000) or depth weighting `W`. Each step is solved by
projected, Jacobi-preconditioned conjugate gradients on the cells off the bounds, with the dense sensitivity
`G` held once in Float32. β starts from the ratio of the largest eigenvalues of the data and model terms, is
cooled to the target misfit, and is bisected when a cooled step overshoots it.

The INR inversion drops the explicit regularizer and parameterises the model as `ρ(x) = mid + half·tanh(f_θ(ψ(x)))`
over Fourier features `ψ` of the cell centres. Its gradient chains the analytic data-misfit gradient `Jᵀr`
with a Zygote pullback of the network only, so the sensitivity is never traced by AD. Adam is not monotone,
so the run keeps the best-misfit network and, on divergence, falls back to it with half the rate.

Training is scheduled rather than flat. The rate ramps up over the first 2 % of the epoch budget, then
decays by a cosine to `Learning rate floor` of its peak, and the Fourier bands are released one at a time
over the first `Frequency annealing` of the budget, lowest first — the coarse-to-fine schedule of the INR
literature, so the network fits the long wavelengths of the anomaly before it has the means to fit the short
ones. Both are on by default (`Learning rate schedule : cosine`, `Frequency annealing : 0.5`) and `none`/`0`
restores the flat training of earlier releases.

What the schedule buys is a better model at the same misfit, mostly in how much mass is recovered and where
it sits — not a better fit. Over six seeds on `demo2`, holding the network fixed, the correlation with
the truth is unchanged within the seed spread (0.77 → 0.78 ± 0.03) while the recovered mass of the compact
body goes from 2.95 Mt to 2.16 Mt against a true 2.01 Mt and its centroid from 201 m to 136 m against a true
147 m, at about 30 % more epochs. On `demo1`, over three seeds, the correlation is likewise unchanged
and the recovered deficit and excess move from 7.0 and 4.9 kt to 5.9 and 3.0 kt against a true 5.0 and
2.5 kt. The schedule also keeps training stable at a hotter peak rate, so the network can be run wider and
faster — which is where the rest of `demo2`'s margin comes from.

One practical point: the network should parameterise the region of interest, not the mesh padding. Nothing
in the INR pays for the volume of a padding cell — no smallness, no depth weighting — so with the padding
free the network fills those oversized far cells to buy misfit. On `examples/demo1` that stalls the run
at RMS 1.35, while the same network on the core alone reaches the target in seconds and correlates better
with the truth. `Invert padding : no` in the control file holds the padding fixed, and it is the default for
`InvCtrl.INR`; `MakeMesh3D(...; invert_padding = false)` bakes the same restriction into the mask instead.

## Status

Forward modelling (prism and FD), model/data/control I/O, synthetics, both inversions and mesh design
(`MakeMesh3D`) are implemented and tested (`test/`); the plotting layer is implemented and exercised by the
two worked synthetics in `examples/`, but has no tests of its own. The field-data workflow is not written
yet.

## References

- Plouff, D., 1976, Gravity and magnetic fields of polygonal prisms and application to magnetic terrain corrections: Geophysics, 41, 727–741.
- Blakely, R. J., 1995, Potential Theory in Gravity and Magnetic Applications: Cambridge University Press.
- Li, Y., and Oldenburg, D. W., 1996, 3-D inversion of magnetic data: Geophysics, 61, 394–408.
- Li, Y., and Oldenburg, D. W., 1998, 3-D inversion of gravity data: Geophysics, 63, 109–119.
- Li, Y., and Oldenburg, D. W., 2000, Joint inversion of surface and three-component borehole magnetic data: Geophysics, 65, 540–552.
- Nagy, D., Papp, G., and Benedek, J., 2000, The gravitational potential and its derivatives for the prism: Journal of Geodesy, 74, 552–560.

## License

MIT License — see [LICENSE.md](LICENSE.md).
