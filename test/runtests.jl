# Author: @pankajkmishra

using Test
using GravityGeophysics

@testset "GravityGeophysics.jl" begin
    
    @testset "Mesh Creation" begin
        # Test uniform mesh
        mesh = create_mesh(nx=10, ny=10, nz=5, lx=1000.0, ly=1000.0, lz=500.0)
        @test mesh.nx == 10
        @test mesh.ny == 10
        @test mesh.nz == 5
        @test mesh.lx ≈ 1000.0
        @test length(mesh.dx) == 10
        @test all(mesh.dx .≈ 100.0)
        
        # Test cell centers
        xc, yc, zc = get_cell_centers(mesh)
        @test length(xc) == 10
        @test xc[1] ≈ -450.0  # Centered mesh
        @test zc[1] ≈ 50.0    # First z center
        
        # Test cell edges
        xe, ye, ze = get_cell_edges(mesh)
        @test length(xe) == 11
        @test xe[1] ≈ -500.0
        @test xe[end] ≈ 500.0
    end
    
    @testset "Voxel And XYZ I/O" begin
        test_dir = mktempdir()

        mesh = create_mesh(nx=10, ny=10, nz=5, lx=1000.0, ly=1000.0, lz=500.0)
        model = create_model(mesh)
        add_block_anomaly!(model;
            xmin=-200.0, xmax=200.0,
            ymin=-200.0, ymax=200.0,
            zmin=100.0, zmax=300.0,
            drho=500.0)

        receivers = create_receivers(mesh; nrx=5, nry=5, z=-10.0)
        data = forward(model, receivers, Prism())

        vox_file = joinpath(test_dir, "test_model.vox")
        save_model_vox(vox_file, model)
        @test isfile(vox_file)
        model_loaded = load_model_vox(vox_file)
        @test size(model_loaded.density) == size(model.density)
        @test maximum(abs.(model_loaded.density .- model.density)) < 1e-10

        xyz_file = joinpath(test_dir, "test_data.xyz")
        save_data_xyz(xyz_file, data)
        @test isfile(xyz_file)
        data_loaded = load_data_xyz(xyz_file)
        @test data_loaded.receivers.n == data.receivers.n
        @test maximum(abs.(data_loaded.gz .- data.gz)) < 1e-12

        recv_file = joinpath(test_dir, "test_receivers.xyz")
        save_receivers_xyz(recv_file, receivers)
        @test isfile(recv_file)

        rm(test_dir; recursive = true)
    end
    
    @testset "Forward Modeling - FD" begin
        # Small test mesh
        mesh = create_mesh(nx=8, ny=8, nz=4, lx=800.0, ly=800.0, lz=400.0)
        
        # Model with central anomaly
        model = create_model(mesh)
        add_block_anomaly!(model;
            xmin=-100.0, xmax=100.0,
            ymin=-100.0, ymax=100.0,
            zmin=100.0, zmax=200.0,
            drho=1000.0)
        
        # Receivers
        receivers = create_receivers(mesh; nrx=5, nry=5, z=10.0)
        
        # Forward modeling
        data, info = forward_fd(model, receivers)
        
        @test length(data.gz) == 25
        @test all(isfinite.(data.gz))
        @test info["iterations"] > 0
    end
    
    @testset "Method Comparison" begin
        # Larger mesh for meaningful comparison
        mesh = create_mesh(nx=16, ny=16, nz=8, lx=1600.0, ly=1600.0, lz=800.0)
        
        model = create_model(mesh)
        add_block_anomaly!(model;
            xmin=-200.0, xmax=200.0,
            ymin=-200.0, ymax=200.0,
            zmin=200.0, zmax=400.0,
            drho=500.0)
        
        receivers = create_receivers(mesh; nrx=9, nry=9, z=20.0)
        
        # Both methods
        data_prism = forward(model, receivers, Prism())
        data_fd, _ = forward_fd(model, receivers)
        
        # Compute stats
        stats = compute_stats(data_prism, data_fd)
        
        # Results should be reasonably close
        @test stats["rel_rmse_pct"] < 10.0  # Within 10%
    end
    
    @testset "Constants" begin
        @test G_NEWTON ≈ 6.67430e-11
        @test SI_TO_MGAL == 1e5
    end

    @testset "Synthetic Demo Bundle" begin
        outdir = mktempdir()
        bundle = generate_demo_bundle(outdir)
        @test haskey(bundle, :model)
        @test haskey(bundle, :data)
        @test isfile(bundle[:model_vox_path])
        @test isfile(bundle[:data_xyz_path])
        @test !any(endswith(name, ".msh") || endswith(name, ".den") || endswith(name, ".obs") for name in readdir(outdir))

        bench = generate_block_benchmark_bundle(outdir)
        @test haskey(bench, :noise_std_si)
        @test isfile(bench[:model_vox_path])
        @test bench[:receivers].n == bench[:mesh].nx * bench[:mesh].ny

        pyhasalmi_mesh = create_pyhasalmi_mesh(dx = 200.0, dy = 200.0, dz = 120.0, lx = 1200.0, ly = 1200.0, lz = 720.0)
        pyhasalmi = generate_pyhasalmi_bundle(outdir; mesh = pyhasalmi_mesh, noise_std_mgal = 0.01, seed = 7, receiver_kwargs = (nrx = 7, nry = 7, xfrac = 0.7, yfrac = 0.7))
        @test haskey(pyhasalmi, :host_density)
        @test isfile(pyhasalmi[:model_vox_path])
        @test isfile(pyhasalmi[:data_xyz_path])
    end

    @testset "Inversion Utilities" begin
        mesh = create_mesh(nx=6, ny=6, nz=4, lx=600.0, ly=600.0, lz=400.0)
        model = create_model(mesh)
        add_block_anomaly!(model; xmin=-100.0, xmax=100.0, ymin=-100.0, ymax=100.0, zmin=100.0, zmax=250.0, drho=200.0)
        receivers = create_receivers(mesh; nrx=4, nry=4, z=-10.0)
        data = forward(model, receivers, Prism())
        G = sensitivity_matrix_prism(model, receivers)
        @test size(G) == (receivers.n, mesh_cell_count(mesh))

        inverted, info = invert_tikhonov(data, mesh; G = G, lambda = 0.5)
        @test size(inverted.density) == size(model.density)
        @test haskey(info, "rmse_mgal")

        inr_model, inr_info = invert_inr(data, mesh; G = G, epochs = 3, hidden_dim = 8, depth = 2, nfreq = 2, verbose = false)
        @test size(inr_model.density) == size(model.density)
        @test length(inr_info["loss_history"]) == 3
    end
    
end
