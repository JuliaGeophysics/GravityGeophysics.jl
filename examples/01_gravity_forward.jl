# Author: @pankajkmishra
#=
================================================================================
Example: 3D Gravity Forward Modeling with GravityGeophysics.jl
================================================================================

This example demonstrates:
1. Creating a mesh
2. Building a density model with block anomaly
3. Computing gravity using Prism (analytical) and FD (numerical) methods
4. Comparing results and saving in UBC format

Test Case:
- Single cubic density anomaly (Δρ = 300 kg/m³)
- Cube size: 8 km × 8 km × 8 km, centered at depth 12 km
- Domain: 160 km × 160 km × 80 km
- Grid: 64 × 64 × 40 cells

Usage (from package root):
    julia --project=. examples/01_gravity_forward.jl
    
Or from examples folder:
    julia --project=.. 01_gravity_forward.jl
    
================================================================================
=#

# Activate the parent package
using Pkg
script_dir = @__DIR__
Pkg.activate(joinpath(script_dir, ".."))

using GravityGeophysics
using Printf
using Dates

function main()
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    run_dir = joinpath(pwd(), "run_" * timestamp)
    mkpath(run_dir)
    println("Run directory: " * run_dir)

# -----------------Create Mesh-------------------------------------
    println("Creating mesh...")

    mesh = create_mesh(
        nx = 64,
        ny = 64,
        nz = 40,
        lx = 160000.0,  # 160 km
        ly = 160000.0,  # 160 km
        lz = 80000.0    # 80 km depth
    )
    
    mesh_info(mesh)

    mesh_path = joinpath(run_dir, "mesh.msh")
    save_mesh_ubc(mesh_path, mesh)

# -----------------Create Density Model---------------------------
    println("Creating density model...")

    model = create_model(mesh; background=0.0)
    
    # Cubic anomaly centered at (0, 0, 12km), size 8km³
    cube_half = 4000.0
    cube_depth = 12000.0
    drho = 300.0
    
    add_block_anomaly!(model;
        xmin = -cube_half, xmax = cube_half,
        ymin = -cube_half, ymax = cube_half,
        zmin = cube_depth - cube_half, zmax = cube_depth + cube_half,
        drho = drho
    )
    
        # Count anomaly cells
        n_anomaly = count(!iszero, model.density)
        @printf("Anomaly: Δρ=%.1f kg/m³, cells=%d (%.1f%%)\n", drho,
            n_anomaly, 100.0 * n_anomaly / (mesh.nx * mesh.ny * mesh.nz))

        model_path = joinpath(run_dir, "density_model.den")
        save_model_ubc(model_path, model; header="# Density Model [kg/m³]")

# -----------------Create Receivers--------------------------------
        println("Creating receivers...")

    receivers = create_receivers(mesh; nrx=51, nry=51, z=50.0, xfrac=0.9, yfrac=0.9)
    
        @printf("Receivers: %d (%d x %d), z=%.1f m\n",
            receivers.n, receivers.nx, receivers.ny, receivers.z[1])

# -----------------Forward Modeling: Prism-------------------------
        println("Forward modeling: Prism...")

    t0 = time()
    data_prism = forward(model, receivers, Prism())
    t_prism = time() - t0
    
        @printf("Prism time: %.3f s\n", t_prism)
        @printf("Prism gz range: [%.6f, %.6f] mGal\n",
            minimum(data_prism.gz) * SI_TO_MGAL,
            maximum(data_prism.gz) * SI_TO_MGAL)

# -----------------Forward Modeling: FD----------------------------
        println("Forward modeling: FD...")

    fd_options = FDOptions(tol=1e-10, maxiter=5000, use_prism_bc=false)
    
    t0 = time()
    data_fd, fd_info = forward_fd(model, receivers; options=fd_options)
    t_fd = time() - t0
    
    @printf("FD time: %.3f s (CG iters: %d)\n", t_fd, fd_info["iterations"])
    @printf("FD gz range: [%.6f, %.6f] mGal\n",
            minimum(data_fd.gz) * SI_TO_MGAL,
            maximum(data_fd.gz) * SI_TO_MGAL)

# -----------------Compare Results---------------------------------
    println("Comparing methods...")

    stats = compute_stats(data_prism, data_fd)

    stats_path = joinpath(run_dir, "stats.txt")
    open(stats_path, "w") do io
        @printf(io, "RMSE: %.6e m/s^2 (%.6f mGal)\n", stats["rmse_si"], stats["rmse_mgal"])
        @printf(io, "MAE: %.6e m/s^2 (%.6f mGal)\n", stats["mae_si"], stats["mae_mgal"])
        @printf(io, "MAX: %.6e m/s^2 (%.6f mGal)\n", stats["max_si"], stats["max_mgal"])
        @printf(io, "Rel RMSE: %.3f%%\n", stats["rel_rmse_pct"])
        @printf(io, "Rel MAX: %.3f%%\n", stats["rel_max_pct"])
        @printf(io, "Prism time: %.3f s\n", t_prism)
        @printf(io, "FD time: %.3f s\n", t_fd)
        @printf(io, "FD iterations: %d\n", fd_info["iterations"])
    end

# -----------------Save Outputs------------------------------------
    println("Saving outputs...")

    data_prism_path = joinpath(run_dir, "gravity_prism.obs")
    data_fd_path = joinpath(run_dir, "gravity_fd.obs")
    save_data_ubc(data_prism_path, data_prism; 
                  header="! Gravity Data - Prism Method [gz in m/s²]")
    save_data_ubc(data_fd_path, data_fd; 
                  header="! Gravity Data - FD Method [gz in m/s²]")

    plot_path = joinpath(run_dir, "gravity_comparison.png")

# -----------------Create Plot-------------------------------------
    if GravityGeophysics.has_visualization
        try
            plot_gravity_comparison(data_prism, data_fd;
                                   title_prefix="Gravity gz",
                                   save_path=plot_path)
        catch e
            @warn "Could not create plot: $e"
        end
    end

    summary_path = joinpath(run_dir, "summary.txt")
    open(summary_path, "w") do io
        println(io, "Run directory: " * run_dir)
        println(io, "Mesh: " * mesh_path)
        println(io, "Model: " * model_path)
        println(io, "Prism data: " * data_prism_path)
        println(io, "FD data: " * data_fd_path)
        println(io, "Stats: " * stats_path)
        println(io, "Plot: " * plot_path)
    end

    println("Done.")

    return data_prism, data_fd, stats
end

# Run main
data_prism, data_fd, stats = main()
