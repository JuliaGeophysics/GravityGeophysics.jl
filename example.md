# GravityGeophysics Example Workflow

This package now uses two exchange formats only:

- Voxel models: `.vox`
- Gravity data and receiver tables: `.xyz`

## Load A Model And Data File

```julia
using GravityGeophysics

model = load_model_vox("examples/demo_output/inversion_YYYYMMDD_HHMMSS/true_model.vox")
data = load_data_xyz("examples/demo_output/inversion_YYYYMMDD_HHMMSS/pyhasalmi_observed_gravity.xyz")
```

## Open The 3D Scene

```julia
using GravityGeophysics

model = load_model_vox("path/to/model.vox")
data = load_data_xyz("path/to/data.xyz")

launch_model_viewer(
    model;
    data = data,
    shapefile_paths = String[],
    open_fullscreen = false,
    block = true,
)
```

## TerraScope-Like Things You Can Do Here

- View horizontal density slices through the model.
- Toggle receiver display on top of the model.
- Show an isovolume for stronger density contrasts.
- Export a viewer snapshot as a PNG.
- Overlay shapefiles with `shapefile_paths = ["path/to/boundary.shp"]`.

## Make High-Quality Static Plots

```julia
using GravityGeophysics

model = load_model_vox("path/to/model.vox")
data = load_data_xyz("path/to/data.xyz")

plot_gravity_data(data; title = "Observed Gravity", save_path = "observed_gravity.png")
plot_model_slice(model; slice_axis = :z, slice_index = 6, save_path = "density_slice.png")
```

All static figures are written with CairoMakie and use original easting and northing coordinates in meters.

## Generate A Self-Contained Example Dataset

```julia
using GravityGeophysics

bundle = generate_pyhasalmi_bundle("examples/demo_output/my_case"; noise_std_mgal = 0.025, seed = 42)

model = load_model_vox(bundle[:model_vox_path])
data = load_data_xyz(bundle[:data_xyz_path])

launch_model_viewer(model; data = data, block = true)
```

The included Pyhasalmi-style case is centered near `E = 452528.63 m`, `N = 7059329.16 m` in TM35FIN and writes only `.vox`, `.xyz`, `.png`, and summary files.

## Run The Included Examples

```julia
julia --project=. examples/01_gravity_forward.jl
julia --project=. examples/02_inversion_comparison.jl
julia --project=. examples/04_pyhasalmi_inversion.jl
julia --project=. examples/03_model_viewer.jl
```

`examples/02_inversion_comparison.jl` is the Python-style dipping-block benchmark with orthogonal `XY`, `XZ`, and `YZ` comparison plots.

`examples/04_pyhasalmi_inversion.jl` is the realistic Pyhasalmi-area synthetic case with the same orthogonal comparison layout plus bulk-density voxel outputs.

The saved example folders contain `.vox`, `.xyz`, `.png`, and summary text files only.