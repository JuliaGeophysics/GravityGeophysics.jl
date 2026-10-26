# Neural inversion (INR)

An **implicit neural representation** replaces the usual "one unknown per cell" with a small
neural network that takes a position and returns the density there:

```
ρ(x, y, z) = network(x, y, z)
```

The unknowns are the network's weights, not the cells. Run it exactly like Gauss–Newton, with
an `InvCtrl.INR` in place of the `InvCtrl.GN`:

```julia
run = Invert3D("model.start", "data.xyz", "FwdCtrl",
               "InvCtrl.INR", "C3.dat", "model.prior")
```

## Why this helps

Two things follow from the change of parameterization.

**The model is a continuous function of position.** Real geology often is too — a basin that
deepens gradually, an alteration halo that fades out, a weathering profile. A cell-based model
has to approximate that with a staircase; a coordinate network represents it directly.

**The regularization is the network, not a penalty term.** There is no smallness or smoothing
term here at all. A small network simply *cannot* produce a speckled model — it does not have
the capacity. Instead of tuning weights on a penalty, you choose how much capacity to give it.

The density is squashed through a `tanh` into your `Density bounds`, so bounds are respected
by construction rather than by projection.

## The settings that matter

```
Hidden width        : 128
Hidden layers       : 3
Fourier frequencies : 4
Learning rate       : 0.015
Frequency annealing : 0.5
Maximum number of iterations : 1500
```

**`Maximum number of iterations` counts training epochs**, so it wants hundreds to thousands,
not tens. This is the biggest difference from Gauss–Newton in day-to-day use.

**`Hidden width` and `Hidden layers`** are the capacity, and therefore the regularization.
Start around 128 × 3. Too small and the model cannot represent the structure; too large and it
starts fitting noise.

**`Fourier frequencies`** set how fine a spatial feature the network can express. A plain
coordinate network is biased towards smooth functions and struggles with sharp edges; Fourier
features of the coordinates fix that. More frequencies means sharper structure is available.

**`Frequency annealing`** releases those frequency bands coarse-to-fine over this fraction of
the training budget. The network finds the broad shape first and adds detail later, which is
both more stable and closer to how you would interpret the data yourself.

**`Learning rate schedule`** ramps up over the first 2 % of the budget, then decays by a cosine
to `Learning rate floor` of the peak. This matters more than it sounds: with the schedule on,
the same network reaches a better model and tolerates a hotter peak rate.

**`Seed`** fixes the network initialisation. Change it and re-run as a sanity check — if the
answer moves a lot, you are asking for more detail than the data support.

## Keep the padding out

```
Invert padding : no
```

This is the default for `InvCtrl.INR` and you should leave it alone.

Nothing in the INR pays for the *volume* of a cell — there is no smallness, no depth
weighting. The padding cells around your survey are enormous compared to core cells, so if
they are free the network discovers it can buy misfit very cheaply by filling them with
far-field mass, and the run stalls well short of the target. Holding the padding fixed, the
same network reaches the target quickly and correlates better with the truth.

## What training looks like

Adam is not monotone — the misfit goes down overall but wanders on the way. Two things handle
that:

- the **best** network seen is always kept, so a late bad epoch cannot cost you the result;
- if the misfit blows past `Divergence factor` times the best, training falls back to the best
  network and halves the learning rate, up to `Maximum restarts` times.

The gradient chains the analytic data-misfit gradient with an automatic-differentiation
pullback of the network only — the sensitivity matrix is never traced by AD, which is what
keeps it affordable.

## Strengths and limits

**Good at**: smoothly varying structure; sparse or irregular surveys, where there is little
data and the smoothness of the representation does the interpolating; getting depth right,
because it has no incentive to pile mass just under the stations.

In `examples/demo2` — a graded basin under 150 scattered stations — it recovers a
−2.89 Mt deficit centred at 54 m against a true −3.05 Mt at 52 m, where Gauss–Newton gives
−6.34 Mt at 189 m.

**Less good at**: sharp, blocky contacts, where the smooth representation rounds corners; and
it is slower in wall-clock terms, needing hundreds of epochs where Gauss–Newton needs a few
iterations. On `examples/demo1`, a blocky truth under a dense regular grid, Gauss–Newton
edges it (correlation 0.66 against 0.63).

See [Choosing between them](choosing.md).
