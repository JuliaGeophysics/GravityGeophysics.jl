# Plots

The plotting layer writes PNG files with [CairoMakie](https://docs.makie.org), so it works
headless — on a server, in CI, over SSH — with no display attached.

## Model viewers

A model opens in an interactive window, one view per entry point:

```julia
PlotModelXY("model.den", "data.xyz")                     # depth slices
PlotModelXZ("model.den", "data.xyz")                     # north-south sections
PlotModelYZ("model.den", "data.xyz")                     # east-west sections
PlotModelXYZ("model.den", "data.xyz")                    # the three planes in 3D
```

`model_file` can also be a list of `label => model` pairs (a path or an already-read
`GravModel3D`), drawn side by side on one slider and one colour scale — the usual way to put
a true model and the models recovered by two algorithms on the same page:

```julia
models = ["True" => "model.true", "GN" => gn.final, "INR" => inr.final]
fig, parts = PlotModelXY(models, "data.xyz"; outline_file = "model.true", interactive = false)
parts.export_figure("ModelXY.png")
```

Each window has a slider with Prev/Next and First/Last, a core/padding toggle, a zoom reset
and an Export Figure button that writes a timestamped PNG into `output_dir`. They return
`(fig, parts)`, so a script can drive the slider and export without touching the window:

```julia
fig, parts = PlotModelXY("model.den", "data.xyz"; interactive = false)
set_close_to!(parts.slider, 12)
parts.export_figure()
```

- **`data_file`** — the second positional argument, drawn as the stations, which immediately
  shows you which structure is actually under coverage; it also sets the axis labels from the
  survey's CRS.
- **`outline_file`** — a model on the same mesh; its nonzero cells are outlined on every slice.
- **`with_padding`** — `false` by default, so you see only the core. The padding is a
  numerical device, not a result; showing it makes the survey look bigger than it was.
- **`max_depth`** — cut the plot off below the depth you trust.
- **`density_range`** — fix the colour scale, so two runs can be compared directly.
- **`show_grid`** — draw the cell edges over the slice.
- **`interactive`** — `false` builds the figure and returns without opening a window.

Air is left blank and sections carry the topography as a staircase, so you can see where the
ground actually is.

## Data maps

```julia
PlotData3D("data.xyz")                                   # observed only
PlotData3D("data.xyz"; predicted_file = "data.pred")     # observed, predicted, residual
```

The residual panel is in **mGal**, on its own symmetric colour scale, with the largest
absolute residual in the title. A residual that is a few per cent of the anomaly and looks
like noise is what you want; structure left in it means the model has not explained
something.

## Convergence

```julia
plot_inv3d_convergence(run.history; output_path = "convergence.png", target_rms = 1.0)
```

RMS against iteration, with the regularization weight β beside it when the algorithm cools
one — β is zero throughout an INR run, since there is no explicit regularizer.

Read it for two things: did the curve flatten (if it is still falling steeply, give it more
iterations), and did it reach the target line.

## From the command line

The `examples/` folder has ready-made scripts with a settings block at the top:

```bash
julia --project=. examples/plot_data3D.jl examples/demo1/data.xyz
julia --project=. examples/plot_model3D.jl examples/demo1/model.true xy
```

## Colours

The default is a reversed `:Spectral`, diverging about zero, which suits density contrast —
positive one way, negative the other, background pale. Pass `colormap` to change it:

```julia
using CairoMakie: Reverse
PlotModelXY("model.den"; colormap = Reverse(:RdBu))
```

Keep a diverging map centred on zero for contrast models. A sequential map hides the sign,
which is the most important thing in a density-contrast figure.
