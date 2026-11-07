# Building a mesh

You rarely want to design a gravity mesh by hand. `MakeMesh3D` reads your survey and writes
everything the inversion needs to start:

```julia
MakeMesh3D("data.xyz")
```

That writes two files next to the data:

- `mesh_start_model.den` — the mesh, filled with your background density contrast. Use it as
  the start model, and as the prior too.
- `C3.dat` — the covariance, with air cells fixed and the rest free.

## What it decides for you

A gravity mesh has three parts, and `MakeMesh3D` sizes each from the survey itself:

**The core** sits under your stations with uniform cells. The cell size comes from how far
apart your stations are — the larger of the median nearest-neighbour distance and
`sqrt(area / n)`. There is no point in cells much smaller than your station spacing; gravity
simply cannot see that detail.

**The padding** grows geometrically outwards. It is there so the edge of the mesh does not
pretend the earth stops at the edge of your survey. Padding cells are big — tens of times the
volume of a core cell — which matters later when you decide whether to invert them.

**The layers** start thin at the surface and grow with depth, down to a depth set as a
fraction of the narrower core width. Gravity resolution falls off with depth, so fine layers
deep down would only be buying you nothing.

The printout tells you what you got:

```
289 data, spacing 94 m mean, 100 m nearest; sensors 4 to 33 m; ground 4 to 38 m
grid 44 x 44 x 23 (44528 cells, 3234 air); core 32 x 32 @ 50 m; first layer 10 m;
  pad 1.14 km/side; depth 0.89 km; G ~ 0.04 GB
```

That last number, `G ~`, is how much memory the sensitivity matrix will take. It is built in memory, so if it is larger
than your machine likes, make the cells coarser or the mesh shallower.

## The knobs

```julia
MakeMesh3D("data.xyz";
           out_model       = "model.start",
           out_covariance  = "C3.dat",
           topography      = (north, east) -> elevation_in_metres,
           cell_width_frac = 1.0,    # core cell as a fraction of station spacing
           n_pad           = 6,      # padding cells per side
           pad_factor      = 1.4,    # how fast they grow
           first_layer_div = 2.0,    # top layer = core cell / this
           vertical_factor = 1.15,   # layer growth with depth
           depth_mult      = 0.5,    # model depth as a multiple of core width
           background_density = 0.0)
```

The two you will reach for most often:

- **`cell_width_frac`** — drop it to `0.5` for cells half the station spacing when you want
  more detail and can afford the memory.
- **`first_layer_div`** — raise it when you have real topography. The relief has to be
  resolved by the top few layers; if your top layer is thicker than your relief, the mesh
  cannot represent the terrain and `MakeMesh3D` will carve the ground down to a layer boundary
  to keep stations in air.

## Topography

Pass a function of position and cells above the ground become air:

```julia
MakeMesh3D("data.xyz"; topography = (north, east) -> dem_lookup(north, east))
```

Without it the ground is flat, at the lowest station elevation by default or wherever
`surface_elevation` says.

Every station has to end up in an air cell (or above the mesh entirely). `MakeMesh3D` checks
this and lowers the ground in any column where a station would otherwise be buried, reporting
each adjustment. A lot of adjustments is a sign your top layers are too thick for your
terrain.

## Finding the core later

The padding is a mesh detail, not something you want in a figure or a mass total. These
helpers find the core from the mesh itself, so they cannot drift out of step with it:

```julia
ix, iy = lateral_core_ranges(m)   # index ranges of the core in x and y
core = m.A[ix, iy, :]             # core only
```

The plotting and scoring routines use the same ranges, which is why a model and its figures
always agree about where the survey was.
