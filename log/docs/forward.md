# Forward modelling

Forward modelling answers one question: **if the ground looked like this, what would the
gravimeter read?** Everything else in the package is built on it.

```julia
using GravityGeophysics

model = ReadModel3D("model.true")
data  = ReadData3D("data.xyz")

pred = forward3d(model, data)     # predicted gz in mGal, one per station
```

`pred` lines up with `data.value` station by station, so comparing them is immediate.

## From files

The file equivalent writes a `.pred` file with the same station positions and errors as the
input:

```julia
ForwardSolve3D("model.true", "data.xyz", "FwdCtrl")
```

## How it works

Each cell is a rectangular prism of uniform density contrast, and the vertical attraction of
such a prism has a closed form (Plouff 1976; Blakely 1995; Nagy et al. 2000). The solver
evaluates it once per mesh node and differences over the eight corners of every cell, which is
exact — there is no discretisation error in the gravity calculation itself, only in how well
your mesh represents the real earth.

Following Blakely, `γ` is the gravitational constant and `ρ` the density contrast. Air cells
(`NaN`) carry no mass. `gz` is positive downwards, so a dense body gives a positive anomaly
when it is below the instrument and a negative one when it is above.

The work is threaded over the data, so start Julia with `-t auto`.

## Where the stations have to be

A station must be in air — either in a cell the model marks as air, or above the top of the
mesh entirely. If one ends up inside rock, you get an error naming the station:

```
146 observation(s) lie inside earth cells, first datum 6 at x = 6672682.2, y = 379397.8,
elevation 27.0 m (depth 11.1 m below the mesh top); raise the sensors or rebuild the mesh
with the topography
```

Two usual causes: you built the mesh without topography when your survey has relief, or your
top layers are too thick to resolve the terrain. See [Building a mesh](data/mesh.md).

```julia
check_receivers(model, data)   # returns nothing if every station is fine
```

## The sensitivity matrix

Inversion needs the derivative of every datum with respect to every cell. That matrix is the
expensive object in a gravity problem — `n_data × n_cells` — and it is built once and reused:

```julia
G = Sensitivity3D(model, data, FwdCtrl3D())
```

It is held transposed in Float32 so each datum's row is contiguous. `MakeMesh3D` prints the
expected size; if it will not fit comfortably in memory, coarsen the mesh. With
`Write Frechet derivative : yes` it is saved as `sensitivity.jld2` and read back with
`LoadSensitivity3D`.

## The finite-difference cross-check

A second, completely independent solver is included. Rather than summing prisms, it solves
Poisson's equation `∇²Φ = 4πγρ` on the mesh and differentiates the potential:

```julia
fd = forward3d_fd(model, data)
```

It is there to catch mistakes: if the two solvers agree, the physics is almost certainly
right. They agree to about 1–2 % RMS on uniform meshes and under 10 % on strongly stretched
ones, the difference being discretisation in the Poisson solve, not error in the prism sum.

It takes any tensor mesh — stretched, laterally padded, with topography — because the model's
grid is padded with air on every face so that the stations, which sit *above* the mesh rather
than inside it, fall within the solved region.

```julia
fd = forward3d_fd(model, data, FDCtrl3D(tol = 1e-11, maxiter = 20_000, n_pad = 8))
```

For production work use `forward3d`; it is faster and exact.

## Building a synthetic

Useful for testing an idea before you commit a survey to it:

```julia
m = ReadModel3D("model.start")
truth = model_like(m, copy(m.A))

add_prism!(truth, (north1, north2), (east1, east2), (depth1, depth2), 250.0)
add_dyke!(truth; top_centre = (n, e, 50.0), strike = 30, dip = 70,
          length = 800, width = 60, depth_extent = 400, density_contrast = 300.0)

clean = data_like(data, forward3d(truth, data))
noisy = add_noise(clean; floor = 0.005, percent = 2.0)   # writes the noise as the error
```
