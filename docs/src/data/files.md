# Survey and model files

Everything the package reads and writes is plain text, so you can open any of it in an editor
and see what is going on. There are three kinds of file: your **survey**, a **model**, and the
small **control files** that say how to run.

## The frame

One convention runs through all of it:

- **x is north, y is east**
- **z is depth, positive downwards**, measured from the top of the mesh
- the top of the mesh is the **highest ground** in the model, and its elevation is written in
  the model file header as `# Top elevation (m)`

So a station's depth is `top_elevation - station_elevation`. A station standing on a hilltop
has depth near zero; one in a valley has a positive depth, and it must sit in a cell the model
marks as air.

## Survey files (`.xyz`)

A few `#` header lines, then one row per station, the layout following the same east-first,
error-then-geographic column order as a ModEM data file:

```
# SYNTHETIC gravity survey, Otaniemi, Espoo (not real data)
# CRS : EPSG:3067
# Geographic CRS : EPSG:4326
# easting(m) northing(m) elevation(m) value(mGal) error(mGal) latitude(deg) longitude(deg)
    378897.836    6672682.245     23.173     -0.0370226   0.00562441    60.17302296    24.81716853
    378997.836    6672682.245     23.532     -0.0469962   0.00583423    60.17305262    24.81896948
```

- **easting, northing** in metres. Note the file writes *easting first*, while the code holds
  `x` (north) and `y` (east) — `ReadData3D` sorts that out for you.
- **elevation** in metres above the datum of your CRS.
- **value** is the vertical gravity anomaly `gz` in **mGal**.
- **error** is the uncertainty on that value in mGal. If both the error and the geographic
  columns are missing, every error becomes 1 mGal and you get a warning — in that case set
  `Error floor (mGal)` in the inversion control file, or the inversion has no idea what
  "fitting the data" means.
- **latitude, longitude** in degrees, WGS84. A file that does not carry them gets them derived
  from the easting/northing and the CRS on read, so every survey is georeferenced on its own;
  `has_geographic(d)` says whether a loaded survey has them.

`# CRS : ...` names the coordinate system of the easting/northing columns. A geographic CRS
such as `EPSG:4326` (longitude first) is projected to its UTM zone automatically on reading; a
projected CRS is used as it is. `crs =` on `ReadData3D` overrides the header.

```julia
d = ReadData3D("data.xyz")
d.x, d.y, d.elevation, d.value, d.error   # vectors, all the same length
d.lat, d.lon                               # latitude, longitude (degrees), WGS84
ndata(d)                                   # how many stations
```

## Model files (`.den`)

A model is a **density contrast** in kg/m³ on a rectangular tensor mesh — the difference from
your background rock, not an absolute density. Zero means "ordinary background".

```
# Written by GravityGeophysics.jl WriteModel3D, density contrast (kg/m^3), ...
# Top elevation (m) : 38.109
44 44 23 0 LINEAR
<44 cell widths in x>
<44 cell widths in y>
<23 layer thicknesses in z>
<values, one block per layer>
<origin: north east depth>
```

Cells above the ground are **air**. They are written as `1e17` in the file and come back as
`NaN` in memory, which is how the rest of the package knows to leave them alone — air carries
no mass and is never inverted.

```julia
m = ReadModel3D("model.start")
m.A                    # the 3D array of density contrast, NaN in air
size(m.A)              # (nx, ny, nz)
m.cx, m.cy, m.cz       # cell centres (north, east, depth)
m.dx, m.dy, m.dz       # cell sizes
m.top_elevation        # elevation of depth zero
air_mask(m)            # true where the cell is air
cell_volumes(m)        # m³ per cell, handy for totalling mass
```

Three models usually appear together, all on the same mesh:

| File | Role |
|---|---|
| `model.start` | where the inversion begins, normally all zeros |
| `model.prior` | what the model is pulled back towards; normally the same as the start |
| `model.true` | only in synthetic tests — the answer, for scoring |

## Covariance file (`C3.dat`)

The model covariance. The inversion reads **only its mask**: `0` means the cell is
fixed (air, or anything you want to hold still), anything else means the cell is free to
change. `MakeMesh3D` writes one for you.

Whether the *padding* around your survey is free is not decided here — it is a property of the
algorithm, and lives in the inversion control file as `Invert padding`. See
[Control files](control_files.md).

## Writing files from Julia

```julia
WriteData3D("out.xyz", d)
WriteModel3D("out.den", m; title = "my model")

# build a new model on the same mesh as an old one
m2 = model_like(m, zeros(size(m.A)))

# build new data on the same stations
d2 = data_like(d, predicted_values)
```
