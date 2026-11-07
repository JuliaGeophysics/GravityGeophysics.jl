# Gauss–Newton inversion

The classical approach, and the one to reach for first. It is fast, predictable, and its
behaviour is well understood after thirty years of use on gravity and magnetic data.

```julia
run = Invert3D("model.start", "data.xyz", "FwdCtrl",
               "InvCtrl.GN", "C3.dat", "model.prior")
```

## The idea

Gravity data never determine a density model on their own — infinitely many models fit any
survey equally well. So instead of asking "which model fits the data", we ask for the model
that fits the data **and** is as simple as possible, and we trade the two off:

```
minimise   ½‖Wd(Gm − d)‖²  +  ½β‖R W (m − mref)‖²
```

The first term is the data misfit. The second is the regularizer: how far the model has
strayed from the prior `mref`, and how rough it is. `β` sets the balance. Density bounds are
enforced throughout by projection, so the answer always stays geologically plausible.

## The settings that matter

### Bounds

```
Density bounds (kg/m3) : -400 400
```

Do not leave these wide. They are the cheapest and most physical regularization you have —
set them from the rock types you actually expect, and the inversion cannot invent densities
that do not exist in your area.

### Smallness and smoothing

```
Smallness weight   : 0.04
Smoothing weight x : 1
Smoothing weight y : 1
Smoothing weight z : 1
```

Smallness keeps the model near the prior; smoothing ties neighbouring cells together. The
weights are relative to a core cell, so `smallness = 1/L²` corresponds roughly to a smoothing
length of `L` core cells.

Raise smoothing if the model is speckled. Raise `z` on its own if you are getting vertical
streaks, which gravity is prone to.

### Weighting — do not skip this

```
Sensitivity weighting    : sensitivity
Depth weighting exponent : 2
Depth weighting z0 (m)   : auto
```

Gravity falls off with distance, so a cell just under a station influences the data far more
than a deep one. Left alone, the inversion exploits this and plasters everything into the top
layer, because that is the cheapest way to fit the data.

- **`sensitivity`** weighting (Li & Oldenburg 2000) uses the sensitivity matrix itself to
  compensate. This is the default and usually the right choice.
- **`depth`** weighting (Li & Oldenburg 1998) uses an analytic `(z + z₀)^(-β/2)` instead.
- **`none`** is available but expect a surface-hugging model.

### Beta cooling

```
Initial damping factor lambda : auto
Lambda cooling factor         : 2
```

`β` starts large (a smooth, heavily regularized model) and is divided by the cooling factor
after each accepted iteration, letting structure in gradually until the target misfit is
reached. `auto` picks the starting value from the ratio of the largest eigenvalues of the data
and model terms, which is almost always better than guessing.

If a cooled step overshoots the target misfit, the solver bisects `β` geometrically to land
in `[0.9, 1] ×` target rather than blowing past it — the discrepancy principle.

## Under the hood

Each iteration solves for a step with projected, Jacobi-preconditioned conjugate gradients
over the cells that are off their bounds, then does a backtracking line search. The dense
sensitivity matrix is held once in Float32. `Maximum CG iterations` and `CG tolerance` control
the inner solve; the defaults are rarely worth changing.

## Strengths and limits

**Good at**: blocky, compact targets. Dense, regular surveys. Being predictable — the same
settings give the same answer, and the parameters mean something you can reason about.

**Less good at**: smoothly varying structure, which the piecewise-constant cell
parameterization represents awkwardly; and sparse, irregular surveys, where it tends to put
mass where the stations happen to be. In `examples/demo2` — a graded basin under a
scattered survey — it recovers twice the mass it needs and buries it more than three times too
deep, and that holds across a wide sweep of its regularization settings.

That is the case the [neural inversion](inr.md) was added for.
