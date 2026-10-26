# API

Everything exported by `GravityGeophysics`, grouped by what you would use it for. In the REPL,
`?` in front of any name gives the same documentation.

## Surveys

```@docs
GravData3D
ReadData3D
WriteData3D
data_like
ndata
has_geographic
PlotData3D
```

## Models

```@docs
GravModel3D
ReadModel3D
WriteModel3D
model_like
same_mesh
air_mask
cell_volumes
surface_depth
sensor_depth
```

## Mesh design

```@docs
MakeMesh3D
design_mesh
build_mesh_bundle
data_spacing
```

## Finding the core

```@docs
edges_from_centers
core_indices
z_indices_for_max_depth
lateral_core_ranges
```

## Control files

```@docs
FwdCtrl3D
ReadFwdCtrl3D
WriteFwdCtrl3D
InvCtrl3D
ReadInvCtrl3D
WriteInvCtrl3D
InvCtrlINR3D
ReadInvCtrlINR3D
WriteInvCtrlINR3D
ReadInvCtrl
Cov3D
ReadCov3D
WriteCov3D
```

## Forward modelling

```@docs
forward3d
ForwardSolve3D
GravKernel3D
sensitivity_row
Sensitivity3D
check_receivers
FDCtrl3D
forward3d_fd
```

## Synthetics

```@docs
sphere_anomaly
add_prism!
add_dyke!
add_noise
chi2_rms3d
```

## Inversion

```@docs
Invert3D
Inv3DOptions
Inv3DResult
AbstractInversion3D
model_weights
GaussNewton3D
GaussNewton3DConfig
GaussNewton3DResult
INR3D
INR3DConfig
INR3DResult
```

## Plotting

```@docs
PlotModelXY
PlotModelXZ
PlotModelYZ
PlotModelXYZ
depth_slice_viewer
xz_slice_viewer
yz_slice_viewer
grav_3d_viewer
draw_north_and_scale!
plot_grav3d_data
plot_inv3d_convergence
compute_colorrange
prepare_model_arrays
model_slice
```

## Constants

```@docs
GAMMA
SI_TO_MGAL
```
