<h1 align="center">GravityGeophysics.jl</h1>

<p align="center"><em>A software repository for gravity research and applications.</em></p>

<p align="center">
	<a href="https://github.com/JuliaGeophysics/GravityGeophysics.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaGeophysics/GravityGeophysics.jl/actions/workflows/CI.yml/badge.svg" alt="CI"></a>
	<a href="https://github.com/JuliaGeophysics/GravityGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

3D gravity forward modelling and LSQR inversion on the tensor meshes of MTGeophysics.jl, with topography and lakes.
Models hold the density contrast (kg/m³) against a reduction density, so the data are the (complete) Bouguer anomaly
at that density: the air cells of the mesh carry no mass, and lake water carries 1000 − ρ_ref.
This is an early version.

## Installation

> [!TIP]
> **New to Julia?** Install it first from [julialang.org/install](https://julialang.org/install/), which sets up Julia and its version manager, [juliaup](https://github.com/JuliaLang/juliaup), in one step. Then open a terminal, type `julia`, and follow the steps below.

Install GravityGeophysics.jl into its own environment (press `]` to enter the package manager):

```julia
julia> ]
pkg> activate @gravitygeophysics
pkg> add https://github.com/JuliaGeophysics/GravityGeophysics.jl
```

Then load it in a Julia session or script:

```julia
julia> using GravityGeophysics
```

In a new session, run `pkg> activate @gravitygeophysics` again before `using GravityGeophysics`.

## Examples

The synthetic benchmark: ground stations over two blocks at the Koillismaa Deep Hole, Finland.

```
julia --project=. helpers/benchmarks_3D.jl
julia --project=. --threads=8 examples/run_inv3D.jl
julia --project=. helpers/plot_models_3D.jl examples/synthetic 5000
```

The KDISynthetic example: a synthetic ground survey over the Koillismaa Deep Intrusion, Finland.

```
julia --project=. --threads=8 examples/KDISynthetic/prepare.jl
julia --project=. --threads=8 examples/KDISynthetic/run_KDISynthetic.jl
julia --project=. examples/KDISynthetic/plot_KDISynthetic.jl
```

Contributions are welcome, see [CONTRIBUTING.md](CONTRIBUTING.md).

> [!NOTE]
> **Also check out** other JuliaGeophysics packages: [MTGeophysics.jl](https://github.com/JuliaGeophysics/MTGeophysics.jl) for magnetotellurics, [MagneticGeophysics.jl](https://github.com/JuliaGeophysics/MagneticGeophysics.jl) for magnetics, [Timekeepers.jl](https://github.com/JuliaGeophysics/Timekeepers.jl) for time-series analysis in electromagnetic geophysics, and [TerraScope.jl](https://github.com/JuliaGeophysics/TerraScope.jl) for looking at diverse geoscientific data together.
