```@raw html
---
layout: home

hero:
  name: "GravityGeophysics.jl"
  tagline: A software repository for gravity geophysics research and applications.
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      text: Run an inversion
      link: /inversion/running
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaGeophysics/GravityGeophysics.jl

features:
  - title: Forward modelling
    details: Closed-form gravity of rectangular prisms over any tensor mesh, with topography. A finite-difference Poisson solver comes along as an independent check.
    link: /forward
  - title: Two inversions, one workflow
    details: Gauss-Newton and a neural (INR) inversion share the same six input files, the same driver and the same result. Swap one control file to change method.
    link: /inversion/running
  - title: Meshes from your survey
    details: Point MakeMesh3D at a data file and it designs a padded mesh, a start model and a covariance around your stations.
    link: /data/mesh
  - title: Plots you can hand over
    details: Depth slices, cross-sections, observed-versus-predicted maps and convergence curves, written straight into the run folder.
    link: /plots
  - title: Worked examples
    details: Two synthetic surveys you can run end to end, with the true model shipped so you can score the result yourself.
    link: /examples
---
```

## What is GravityGeophysics.jl?

GravityGeophysics.jl does 3D forward modelling and inversion of gravity data for **density
contrast**. It is part of the [JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics)
and is driven entirely by small plain-text files: a `.xyz` survey, block-layout `.den`
models, `Key : value` control files, and a run folder written next to your data.

If you have a gravity survey and want a density model under it, you are in the right place.

### Two ways to invert

Both read the same files and return the same kind of result, so you can try one and then the
other without changing anything else:

- **Gauss–Newton** — the classical approach. Smallness and smoothing hold the model together,
  and sensitivity or depth weighting stops everything piling up at the surface.
- **Implicit neural representation (INR)** — a small coordinate network
  ([Lux.jl](https://github.com/LuxDL/Lux.jl)) outputs the density of every cell. There is no
  explicit regularizer; the network's limited capacity is what keeps the model smooth.

[Choosing between them](inversion/choosing.md) explains when each one tends to win.

## Install

```julia
pkg> add GravityGeophysics
```

You need Julia 1.13.0. Then head to [Getting Started](getting_started.md) to run your first
forward model and inversion.

## Where to go next

| If you want to… | Read |
|---|---|
| see it work in five minutes | [Getting Started](getting_started.md) |
| understand the file formats | [Survey and model files](data/files.md) |
| build a mesh around your stations | [Building a mesh](data/mesh.md) |
| predict data from a model | [Forward modelling](forward.md) |
| invert real data | [Running an inversion](inversion/running.md) |
| make figures | [Plots](plots.md) |
| see complete worked cases | [Worked examples](examples.md) |
