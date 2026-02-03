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
    
    @testset "Model Creation" begin
        mesh = create_mesh(nx=10, ny=10, nz=5, lx=1000.0, ly=1000.0, lz=500.0)
        
        # Test empty model
        model = create_model(mesh)
        @test size(model.density) == (10, 10, 5)
        @test all(model.density .== 0.0)
        
        # Test with background
        model2 = create_model(mesh; background=2670.0)
        @test all(model2.density .≈ 2670.0)
        
        # Test block anomaly
        model3 = create_model(mesh)
        add_block_anomaly!(model3;
            xmin=-100.0, xmax=100.0,
            ymin=-100.0, ymax=100.0,
            zmin=100.0, zmax=300.0,
            drho=500.0)
        @test any(model3.density .≈ 500.0)
    end
    
    @testset "Receivers" begin
        mesh = create_mesh(nx=10, ny=10, nz=5, lx=1000.0, ly=1000.0, lz=500.0)
        
        receivers = create_receivers(mesh; nrx=11, nry=11, z=10.0)
        @test receivers.n == 121
        @test receivers.nx == 11
        @test receivers.ny == 11
        @test all(receivers.z .≈ 10.0)
    end
    
    @testset "Forward Modeling - Prism" begin
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
        data = forward(model, receivers, Prism())
        
        @test length(data.gz) == 25
        @test all(isfinite.(data.gz))
        # Gravity should be positive (mass below)
        @test maximum(data.gz) > 0
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
    
    @testset "UBC Format I/O" begin
        # Create temp directory for test files
        test_dir = mktempdir()
        
        # Test mesh I/O
        mesh = create_mesh(nx=5, ny=5, nz=3, lx=500.0, ly=500.0, lz=300.0)
        mesh_file = joinpath(test_dir, "test_mesh.msh")
        save_mesh_ubc(mesh_file, mesh)
        @test isfile(mesh_file)
        
        mesh_loaded = load_mesh_ubc(mesh_file)
        @test mesh_loaded.nx == mesh.nx
        @test mesh_loaded.lx ≈ mesh.lx
        
        # Test model I/O
        model = create_model(mesh)
        add_block_anomaly!(model;
            xmin=-50.0, xmax=50.0,
            ymin=-50.0, ymax=50.0,
            zmin=50.0, zmax=150.0,
            drho=100.0)
        
        model_file = joinpath(test_dir, "test_model.den")
        save_model_ubc(model_file, model)
        @test isfile(model_file)
        
        model_loaded = load_model_ubc(model_file, mesh)
        @test size(model_loaded.density) == size(model.density)
        
        # Test data I/O
        receivers = create_receivers(mesh; nrx=5, nry=5, z=10.0)
        data = forward(model, receivers, Prism())
        
        data_file = joinpath(test_dir, "test_data.obs")
        save_data_ubc(data_file, data)
        @test isfile(data_file)
        
        data_loaded = load_data_ubc(data_file)
        @test data_loaded.receivers.n == data.receivers.n
        @test length(data_loaded.gz) == length(data.gz)
        
        # Cleanup
        rm(test_dir; recursive=true)
    end
    
    @testset "Constants" begin
        @test G_NEWTON ≈ 6.67430e-11
        @test SI_TO_MGAL == 1e5
    end
    
end
