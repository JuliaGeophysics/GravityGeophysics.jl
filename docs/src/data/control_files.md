# Control files

Control files are plain `Key : value` text. Unknown keys, duplicate keys and missing required
keys are all errors, so a typo is caught immediately rather than silently ignored. Anything
after `#` is a comment, and any key you leave out takes its default.

## `FwdCtrl` — the forward solver

```
Write Frechet derivative : no
```

- **`Write Frechet derivative`** — `yes` keeps the sensitivity matrix after the run, as
  `sensitivity.jld2` in the run folder, if you want it for something else. Read it back with
  `LoadSensitivity3D`.

## `InvCtrl.GN` — Gauss–Newton

```
Algorithm                         : GN
Initial damping factor lambda     : auto
Exit search when rms is less than : 1
Maximum number of iterations      : 20
Lambda cooling factor             : 2
Density bounds (kg/m3)            : -400 400
Max line search steps             : 12
Smallness weight                  : 0.04
Smoothing weight x                : 1
Smoothing weight y                : 1
Smoothing weight z                : 1
Sensitivity weighting             : sensitivity
Depth weighting exponent          : 2
Depth weighting z0 (m)            : auto
Error floor (mGal)                : 0
Error percent                     : 0
GN damping                        : 0.001
Maximum CG iterations             : 50
CG tolerance                      : 0.001
Invert padding                    : yes
```

The ones that change your answer most:

| Key | What it does |
|---|---|
| `Exit search when rms is less than` | when to stop. `1` means "fit the data to the size of the errors". Below 1 you are fitting noise. |
| `Density bounds (kg/m3)` | hard limits on the density contrast. Set them from your rock types; they do a lot of the regularization work. |
| `Smallness weight` | how strongly the model is pulled back to the prior. Raise it for a quieter model. |
| `Smoothing weight x/y/z` | how strongly neighbours are tied together. Raise `z` to discourage vertical streaks. |
| `Sensitivity weighting` | `sensitivity`, `depth` or `none`. Without one, everything piles up just under the stations. |
| `Error floor` / `Error percent` | the error actually used is the largest of the file error, this floor, and this percentage of the datum. Use them when you do not trust the error column. |

`auto` for lambda lets the solver pick a starting regularization weight from the problem
itself; that is usually the right choice.

## `InvCtrl.INR` — neural inversion

```
Algorithm                         : INR
Exit search when rms is less than : 1
Maximum number of iterations      : 1500
Density bounds (kg/m3)            : -400 400
Error floor (mGal)                : 0
Error percent                     : 0
Hidden width                      : 128
Hidden layers                     : 3
Fourier frequencies               : 4
Learning rate                     : 0.015
Coordinate scaling                : per-axis
Seed                              : 42
Divergence factor                 : 3
Maximum restarts                  : 3
Learning rate schedule            : cosine
Learning rate floor               : 0.05
Frequency annealing               : 0.5
Invert padding                    : no
```

**`Maximum number of iterations` means training epochs here**, so it wants thousands, not
tens. The settings that matter:

| Key | What it does |
|---|---|
| `Hidden width`, `Hidden layers` | the network's capacity, which *is* the regularization. Wider and deeper means more detail and more risk of fitting noise. |
| `Fourier frequencies` | how fine a spatial feature the network can express. More frequencies, sharper structure. |
| `Learning rate` | too high and training swings; too low and it never arrives. |
| `Frequency annealing` | releases the frequency bands coarse-to-fine over this fraction of the budget, so the network finds the broad shape before the detail. |
| `Seed` | the network's initialisation. Change it to check your result is not an accident. |

## `Invert padding`

This key appears in both files and decides whether the geometric padding around your survey is
free to change.

- **Gauss–Newton: `yes`.** Its sensitivity or depth weighting already prices the volume of a
  padding cell, so those cells cost what they are worth.
- **The INR: `no`** (the default for `InvCtrl.INR`). There is no term in the INR that pays for
  volume, so with the padding free the network fills those huge far-away cells to buy misfit
  cheaply and the run stalls.

Because this is a property of the *algorithm* and not of the mesh, both runs can read the same
`C3.dat`. `Summary.txt` records which setting was used.

## Reading and writing them from Julia

```julia
fwd  = ReadFwdCtrl3D("FwdCtrl")
ctrl = ReadInvCtrl("InvCtrl.GN")      # picks GN or INR from the Algorithm line

WriteInvCtrl3D("InvCtrl.GN", InvCtrl3D(algorithm = :gn, lambda = nothing,
                                       target_rms = 1.0, max_iter = 20))
```

`WriteInvCtrl3D` and `WriteInvCtrlINR3D` always write every key, so they are a good way to see
the full set of defaults.
