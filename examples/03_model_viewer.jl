using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GravityGeophysics

function main()
    outdir = joinpath(@__DIR__, "demo_output", "viewer_demo")
    bundle = generate_demo_bundle(outdir)
    model = load_model_vox(bundle[:model_vox_path])
    data = load_data_xyz(bundle[:data_xyz_path])
    launch_model_viewer(model; data = data, shapefile_paths = String[], open_fullscreen = false, block = true)
end

main()