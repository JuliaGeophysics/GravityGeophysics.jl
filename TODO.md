# TODO: GravityGeophysics.jl

Status: `[ ]` open, `[~]` in progress, `[x]` done. Plan confirmed with the user 2026-11-07 (supersedes the 2026-11-06 draft).

## Goal
GravityGeophysics.jl, MagneticGeophysics.jl and MTGeophysics.jl identical in features, style and naming down to the
small details, except where the physics makes it impossible, so a joint inversion (`../JointInversion/`) can drive
them through one interface. Magnetics is done (`../MagneticGeophysics.jl/TODO.md`, Phase 11); gravity follows it.
The magnetic package is the working template, `~/github/MTGeophysics.jl` (v0.7.3) the frozen reference.

## What is there (v0.3.0, clean tree)
- Model (`Model.jl`): `GravModel3D`, `ReadModel3D` / `WriteModel3D`, `model_like`, `same_mesh`, `air_mask`,
  `cell_volumes`, `surface_depth`, `sensor_depth`. MT block layout, LINEAR density contrast (kg/m³), air 1e17 → NaN,
  but depth measured from the mesh top with a `# Top elevation (m)` header line; rotation ignored
- Data (`Data.jl`): `GravData3D`, `ReadData3D` / `WriteData3D`, `.xyz` columns easting northing elevation gz error
  (lat lon), `# CRS :` header, geographic CRS projected to UTM; `data_like`, `has_geographic`. No MT-style data file,
  no `Frame`, no sensor rotation
- Forward (`Fwd3D.jl`, `Fwd3D_FD.jl`): analytic prism gz (Plouff/Nagy), `forward3d`, `Sensitivity3D` (Gᵀ Float32),
  `Save/LoadSensitivity3D` (JLD2), `check_receivers`, `ForwardSolve3D(model, data, FwdCtrl)`; finite-volume Poisson
  solver `forward3d_fd` / `FDCtrl3D`. gz only
- Controls (`Control3D.jl`): `FwdCtrl3D` (`Write Frechet derivative`, as MT 2D), `InvCtrl3D` (GN, MT λ keys),
  `InvCtrlINR3D`, `ReadInvCtrl` dispatch, `Cov3D` / `C3.dat` (mask only used)
- Inversion (`Inv3D.jl`, `Inv3D_GN.jl`, `Inv3D_INR.jl`): `Invert3D` with the `inv3d_*` interface,
  `AbstractInversion3D`, `Inv3DOptions`, `Inv3DResult`, `Inv3DProblem`, `GaussNewton3D(Config/Result)`;
  `INR3D` (Lux/Zygote/Optimisers, own loop); files `Invert3D(start, data, FwdCtrl, InvCtrl, C3.dat, prior)` →
  `model.den`, `data.pred`, `History.log`, `Summary.txt`, `inputs/`. No NLCG, no VFSA, no exported objective functions
- Mesh (`MakeMesh3D.jl`): `MakeMesh3D(data_file)`, `design_mesh`, `build_mesh_bundle`, `data_spacing`; writes start
  model + `C3.dat`. No GUI, no GeoTIFF, no inversion mask, no `Mask3D.jl`
- Synthetics (`Synthetic3D.jl`): `sphere_anomaly`, `add_prism!`, `add_dyke!`, `add_noise`, `chi2_rms3d`
- Plotting: `PlotModelXY/XZ/YZ/XYZ`, viewers (`depth_slice_viewer`, …), `plot_grav3d_data`, `plot_inv3d_convergence`,
  `PlotData3D`; `CoreUtils3D`. No `ShapefileOverlay`, `GeoRef3D`, `PlotDataMap`, `save_figure`
- No `EditModel3D`, `MeshToMesh`, `Rotate.jl`
- Examples: `examples/demo1`, `demo2` (synthetic, Otaniemi, EPSG:3067; GN and INR), run scripts with plots inside
- Tests: `TestIO3D`, `TestControl3D`, `TestForward3D`, `TestSensitivity3D`, `TestMesh3D`, `TestInversion3D`,
  `TestWorkflows3D`; runtests not in MT's sectioned style
- Deps: JLD2, Lux, Zygote, Optimisers, NCDatasets, DelimitedFiles beyond magnetics; no ArchGDAL, Shapefile,
  GeoInterface. `LICENSE.md`, `README.md`, Documenter `docs/` (several pages)
- `../JointInversion/run_joint_xgrad.jl` calls an older API of both packages (`GG.load_model`, `GG.GravForward`,
  `GG.InvProblem`, `MG.InvProblem`) and runs against neither package as it is now

## Decisions for the alignment (confirmed 2026-11-07)
- Template: magnetics as it is now, file for file, including its deviations from MT (LSQR only, `rotate_data`,
  `rms`, quadratic-fit `GeoRef3D`). MT (`~/github/MTGeophysics.jl` v0.7.3, frozen) only where magnetics has nothing.
  Anything better in gravity is logged under "Better in gravity"
- No "modem" anywhere. Names as magnetics: `Model`, `load_model` / `write_model`, `load_mesh`, `air_cells`,
  `active_cells`, `same_mesh`, `load_mask` / `write_mask`, `water_mask`, MT's topography/air masks; data `Data`,
  `ndata`, `load_data_grav` / `write_data_grav`, `COMPONENTS = (:gz,)`, `Frame`, `load_frame`, `to_frame`,
  `north_azimuth`, `drape!`, `check_frame`, `rotate_data`
- Version 0.4.0, breaking, no deprecation shims
- Model file exactly the MT/magnetic one (absolute z down through the origin line, no `Top elevation` line, no
  converter: old files are an error). Values: density contrast in kg/m³, LINEAR, air 1e17 → NaN. Run output `model.den`
- Data: native file in the magnetic layout (`Code GG_Lat GG_Lon X Y Z Component Value Error`, `> gz`, `> [mGal]`,
  `> azimuth`, `> origin`, `> n`; no `> F I D` line). Today's `.xyz`, CSV and Geosoft XYZ read through the import
  path, as magnetics. Errors come from the data file only. gz only; gradients later if needed
- Terrain: relief in the mesh as air cells, data used as given (no Bouguer/terrain tools). With air carrying no mass
  the data must be a (complete) Bouguer anomaly at the reference density; README says so
- Forward as magnetics: `GravForward(m, d; solver = :dense | :matrixfree)`, `forward`, `jvec`, `jtvec`,
  `model_vector`, `model_array`, `ForwardSolve3D`; `FwdCtrl` has `Forward solver` only. Removed: `Write Frechet
  derivative`, `Save/LoadSensitivity3D`, JLD2, NCDatasets (after checking nothing else uses them)
- FD Poisson solver (`Fwd3D_FD.jl`): unexported, used only as a check in `test/TestAccuracy3D.jl`
- Inversion: LSQR only, `Inv3D.jl` + `Inv3D_LSQR.jl` as magnetics, `InvCtrl.LSQR` with magnetics' keys
  (λ as a multiple of the Hessian ratio, `Model weighting`, `Depth exponent`, LSQR keys), `Density bounds (kg/m3) :
  -1000 1000`, `Max density step (kg/m3) : Inf`; per-cell bounds (`lower`/`upper` models, drillhole constraint);
  exported `objective`, `gradient`, `data_hessian`, `reg_hessian`, `phi_d`, `phi_m`, `rms`, `predicted`; run inputs
  `Invert3D(model.start, data.dat, FwdCtrl, InvCtrl.LSQR, model.mask, model.prior)` → `model.den`, `data.pred`,
  `0lsqr.log`, `inputs/`. Removed: GN, INR, `InvCtrlINR3D`, `ReadInvCtrl`, `C3.dat`/`Cov3D`, `Invert padding`,
  `Depth weighting z0 (m)`, `Error floor (mGal)`, `Error percent`, Lux/Zygote/Optimisers
- Mesh: magnetics' `MakeMesh3D` + GUI, GeoTIFF DEM (ArchGDAL) or `surface(x, y)`, `uniform_layers`; start model +
  mask. Lakes (`Lakes3D.jl`, `read_lakes`, mask 9): water fixed at 1000 − ρ_ref, ρ_ref a keyword (default 2670 kg/m³)
- Ported tools: `EditModelByLayers`, `EditModelByDrawing`, `MeshToMesh`; plotting `PlotModelXY/XZ/YZ/XYZ`,
  `PlotDataMap`, `plot_data`, `plot_data_fit`, `write_data_gis`, `save_figure`, `CoreUtils3D`, `ShapefileOverlay`,
  `GeoRef3D`; GLMakie in `try … catch` as magnetics
- Removed: `Synthetic3D.jl` (bodies move into `helpers/benchmarks_3D.jl`), demo1/demo2 and their scripts
- Examples as magnetics: run and tool scripts in `examples/`, controls in `examples/ctrl/3D/`, `examples/synthetic`
  from `helpers/benchmarks_3D.jl` (same geometry as magnetics, density contrast), `examples/KDISynthetic` (gitignored
  as in magnetics): same DEM, lake and truth geometry, ground stations on a 500 m grid on the DEM, noise
  0.02 mGal + 2 %. Densities proposed by Claude, confirmed by the user before any run; `helpers/plot_*_3D.jl`
- Tests in magnetics' `runtests.jl` layout and file set (`TestAccuracy3D` incl. the FD check)
- Repo: keep `.github` (CI on Julia 1.13, Documenter, TagBot, dependabot); add CITATION.cff, CODE_OF_CONDUCT.md,
  CONTRIBUTING.md from magnetics; docs reduced to magnetics' skeleton, old pages moved to `log/docs/`
- Out of scope: `../JointInversion` scripts; a later pass further toward MT
- Process: one phase at a time; tests pass after each, pure renames reproduce the baseline exactly; stop and report.
  No backup copy beyond the clean git tree. Nothing committed

## Tasks: magnetics alignment (Phase 12), done 2026-11-08
- [x] 12.0 Baseline (v0.3.0, 8 threads): demo1 max |gz| 0.746856195 mGal, demo2 0.549753615; FD vs analytic 5.3 % /
      6.5 % (relative max); 152 tests pass
- [x] 12.1 Model file and names as magnetics (`Model`, `load_model`, `write_model`, `load_mesh`, absolute z); `Mask3D.jl`,
      inversion mask replaces `C3.dat`. The new forward on the old demo files (origin z = −top elevation) gives max |gz|
      0.746856193 / 0.549753613 mGal: 2e-9 mGal from the baseline, the old 1e-6 m station nudge (now removed)
- [x] 12.2 Data: native file (5 '>' lines, no `> F I D`), `load_data_grav` / `write_data_grav`, `Frame`, `drape!`
      (height ≥ 0), `check_frame` (stations on or above the earth cells), `rotate_data`, CSV / Geosoft XYZ import
- [x] 12.3 Forward: `GravForward` dense / matrix-free (agree to 1e-15), `forward`, `jvec`, `jtvec`, `ForwardSolve3D`,
      `FwdCtrl` `Forward solver`; FD solver in `test/FwdFD3D.jl` (TestAccuracy3D 4.6: 1.4 % relative rms); JLD2,
      NCDatasets, Lux, Zygote, Optimisers, DelimitedFiles out
- [x] 12.4 `Control3D.jl`, `Inv3D.jl`, `Inv3D_LSQR.jl` as magnetics; `Density bounds (kg/m3)`, `Max density step (kg/m3)`;
      per-cell bounds; exported objective functions; GN, INR, `Synthetic3D.jl` removed
- [x] 12.5 `MakeMesh3D` + GUI + GeoTIFF, `Lakes3D.jl` (water 1000 − ρ_ref, keyword `ρ_ref`, default 2670)
- [x] 12.6 `EditModel3D`, `MeshToMesh` (density)
- [x] 12.7 Plotting and copied files (`PlotData3D`, `PlotModel*`, `CoreUtils3D`, `ShapefileOverlay`, `GeoRef3D`), GLMakie
      via try/catch
- [x] 12.8 `helpers/benchmarks_3D.jl` → `examples/synthetic` (441 ground stations, blocks +300 / +200 kg/m³; LSQR rms 0.90
      in 8 iterations, 12 s); run/tool scripts, `examples/ctrl/3D`; `examples/KDISynthetic` (gitignored): 729 ground
      stations on a 500 m grid, 0.02 mGal + 2 %, hole density log ± 50 kg/m³; LSQR rms 0.95 in 5 iterations, 66 s;
      four figures. demo1/demo2 removed
- [x] 12.9 Tests in magnetics layout: 296 pass; README, CITATION, CoC, CONTRIBUTING, docs skeleton (old pages in
      `log/docs/`, old README in `log/README_v0.3.0.md`), CI on 1.13 kept; version 0.4.0

## Deviations from magnetics (physics or units)
- Stations may lie on the earth cells (ground surveys): `drape!` height ≥ 0, `check_frame` z ≤ ground, MakeMesh3D lowers a
  column only when a station is below its earth top; the kernel takes the exact limit (x log x → 0) on mesh planes
- Default `Model weighting : depth`, `Depth exponent : 2` (Li & Oldenburg 1998 for gravity): with sensitivity weighting
  the synthetic block's density went to the surface (peak layer 1, correlation 0.48); depth weighting puts it in the
  block (layer 3, 0.58)
- `GravForward.G` is one matrix (gz only); magnetics keeps one per functional
- Lakes: water is 1000 − ρ_ref (not 0); KDI hole bounds are ± a value, not ×/÷ a factor
- No inducing-field line in the data file; `COMPONENTS = (:gz,)`, units mGal

## Vertical datum (2026-11-08, at the user's request: the highest ground is z = 0)
- [x] `Model.datum`: elevation (m a.s.l.) of z = 0, on file as "# Datum elevation (m) : value" (full precision); a file
      without the line (MT, magnetic, v0.3) is at sea level. `same_mesh` compares it; masks, editors, MeshToMesh and the
      inversion output carry it
- [x] `MakeMesh3D` (`mesh_model`) writes the start model with z = 0 at the top of the highest earth cell and that
      elevation as the datum; the bundle stays at sea level
- [x] Data files keep z = −elevation (identical to magnetics); stations meet the mesh at `d.z + m.datum`: `GravForward`
      (`zs`), `check_frame`, `drape!`, depth weighting, the 3D viewer, MeshToMesh between datums
- [x] Tests: file round trip, forward unchanged by a datum shift (1e-12), MakeMesh z = 0 at the highest ground, MeshToMesh
      across datums: 308 pass. Synthetic (datum 250 m) and KDISynthetic (343.75 m) regenerated: same data and runs
      (rms 0.903 / 0.950)
- [ ] MagneticGeophysics.jl: the same datum (asked, gravity first), else a gravity model read there is off by the datum

## Section placement (2026-11-08)
- `helpers/plot_models_3D.jl`: slice and sections through the centre of mass of the body (cells ≥ half the largest Δρ
  connected to it, of the true model, else the first recovered one); sections reach 1.5 × its base unless max_depth is given
- KDISynthetic fig2: A–A′ through the hole, B–B′ on the KDI funnel's axis 700 m south of it (hole projected, dashed)

## Review fixes (2026-11-08)
- Rotation is kept (MT meshes are often turned to strike; the joint inversion runs gravity on them); the forward was
  already exact under rotation (TestAccuracy3D 4.9)
- [x] MakeMesh3D: the mesh carries the data's azimuth as its rotation (was always 0: MakeMesh's own output was rejected)
- [x] MeshToMesh: target cell centres taken into the source frame (rotation difference, as rotate_data) and the data
      file turned to the target rotation (was: raw frame coordinates, a body 90° off on a turned target, silently)
- [x] Editors: lake water (mask 9) is kept like air (`mask_file`, else `<stem>.mask` / `model.mask` beside the model on
      the same mesh); before, an edit over a lake overwrote 1000 − ρ_ref and mask 9 then held the wrong value fixed
- [x] Editor defaults are densities (100, 300 kg/m³; colour range ±300) instead of the magnetic 0.01 / 0.05
- Tests: 323 pass. Magnetics has the same three code paths (plus its amplitude NaN at κ = 0), not changed here

## Open items
- KDISynthetic densities are Claude's proposals (gneiss 2700, granite 2630, gabbronorite 2950, peridotite 3150, BIF 3400,
  diabase 2950 kg/m³ against 2670): confirm or replace. The +30 host everywhere adds a ~9 mGal offset to the data;
  a host at the reduction density (0 contrast) would give a cleaner Bouguer anomaly
- KDI ground stations: 294 columns lowered up to 37 m (the DEM at the station above the 25 m layer the column snaps to)
- GLMakie windows only checked headless

