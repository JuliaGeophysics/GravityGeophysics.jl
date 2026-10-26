<h1 align="center">GravityGeophysics.jl</h1>

<p align="center"><em>A software repository for gravity geophysics research and applications.</em></p>

<p align="center">
	<a href="https://juliageophysics.github.io/GravityGeophysics.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation"></a>
	<a href="https://github.com/JuliaGeophysics/GravityGeophysics.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaGeophysics/GravityGeophysics.jl/actions/workflows/CI.yml/badge.svg" alt="CI"></a>
	<a href="https://github.com/JuliaGeophysics/GravityGeophysics.jl/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT license"></a>
</p>

3D forward modelling and inversion of gravity data for density contrast. Two inversions —
Gauss–Newton and an implicit neural representation ([Lux.jl](https://github.com/LuxDL/Lux.jl)) —
share one driver and one plain-text file workflow, so you can swap one control file to change
method. Part of the [JuliaGeophysics ecosystem](https://github.com/JuliaGeophysics).

## Requirements

- [Julia](https://julialang.org) 1.13.0 — the only supported version

## Installation

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

To get started, see the [documentation](https://juliageophysics.github.io/GravityGeophysics.jl/dev/).
It covers the file formats, building a mesh, running both inversions, and two worked synthetic
examples.

## License

MIT License — see [LICENSE.md](LICENSE.md).
