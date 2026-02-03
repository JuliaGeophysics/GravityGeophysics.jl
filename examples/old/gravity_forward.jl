# Julia 1.11 
# Author: Pankaj K Mishra 
# Status: Experimental  

using SparseArrays, LinearAlgebra, Printf, Plots

struct GravityModel
    nx::Int
    ny::Int
    nz::Int
    dx::Float64
    dy::Float64
    dz::Float64
    rho::Array{Float64, 3}
end

function initialize_density_model(nx, ny, nz, dx, dy, dz)
    rho = fill(1.0, nx, ny, nz)
    cube_size_in_cells = Int(1_000 / dx)
    cube_start_x = (nx - cube_size_in_cells) ÷ 2 + 1
    cube_end_x = cube_start_x + cube_size_in_cells - 1
    cube_start_y = (ny - cube_size_in_cells) ÷ 2 + 1
    cube_end_y = cube_start_y + cube_size_in_cells - 1
    cube_start_z = (nz - cube_size_in_cells) ÷ 2 + 1
    cube_end_z = cube_start_z + cube_size_in_cells - 1
    rho[cube_start_x:cube_end_x, cube_start_y:cube_end_y, cube_start_z:cube_end_z] .= 5000.0
    return GravityModel(nx, ny, nz, dx, dy, dz, rho)
end

function load_gravity_model_ubc(filename::String)
    open(filename, "r") do file
        readline(file)
        nx, ny, nz = parse.(Int, split(readline(file)))
        dx, dy, dz = parse.(Float64, split(readline(file)))
        rho = Array{Float64, 3}(undef, nx, ny, nz)
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    rho[i, j, k] = parse(Float64, readline(file))
                end
            end
        end
        return GravityModel(nx, ny, nz, dx, dy, dz, rho)
    end
end

function construct_laplacian(nx, ny, nz, dx, dy, dz)
    N = nx * ny * nz
    A = spzeros(N, N)
    index = (i, j, k) -> (i - 1) * ny * nz + (j - 1) * nz + k
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                idx = index(i, j, k)
                A[idx, idx] = -2 * (1/dx^2 + 1/dy^2 + 1/dz^2)
                if i > 1   A[idx, index(i-1, j, k)] = 1 / dx^2 end
                if i < nx  A[idx, index(i+1, j, k)] = 1 / dx^2 end
                if j > 1   A[idx, index(i, j-1, k)] = 1 / dy^2 end
                if j < ny  A[idx, index(i, j+1, k)] = 1 / dy^2 end
                if k > 1   A[idx, index(i, j, k-1)] = 1 / dz^2 end
                if k < nz  A[idx, index(i, j, k+1)] = 1 / dz^2 end
            end
        end
    end
    return A
end

function solve_poisson_direct(gravity_model::GravityModel)
    nx, ny, nz = gravity_model.nx, gravity_model.ny, gravity_model.nz
    dx, dy, dz = gravity_model.dx, gravity_model.dy, gravity_model.dz
    G = 6.67430e-11
    N = nx * ny * nz
    b = zeros(Float64, N)
    index = (i, j, k) -> (i - 1) * ny * nz + (j - 1) * nz + k
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                b[index(i, j, k)] = -4 * π * G * gravity_model.rho[i, j, k]
            end
        end
    end
    A = construct_laplacian(nx, ny, nz, dx, dy, dz)
    φ_flat = A \ b
    return reshape(φ_flat, nx, ny, nz)
end

function compute_gravity(gravity_model::GravityModel, φ::Array{Float64, 3})
    nx, ny, nz = gravity_model.nx, gravity_model.ny, gravity_model.nz
    dx, dy, dz = gravity_model.dx, gravity_model.dy, gravity_model.dz
    gx = zeros(Float64, nx, ny, nz)
    gy = zeros(Float64, nx, ny, nz)
    gz = zeros(Float64, nx, ny, nz)
    for i in 2:nx-1
        for j in 2:ny-1
            for k in 2:nz-1
                gx[i, j, k] = -(φ[i+1, j, k] - φ[i-1, j, k]) / (2 * dx)
                gy[i, j, k] = -(φ[i, j+1, k] - φ[i, j-1, k]) / (2 * dy)
                gz[i, j, k] = -(φ[i, j, k+1] - φ[i, j, k-1]) / (2 * dz)
            end
        end
    end
    return gx, gy, gz
end

function plot_and_save_results(model::GravityModel, gz::Array{Float64,3})
    slice_idx_z = Int(model.nz / 2)
    density_map = model.rho[:, :, slice_idx_z]
    bouguer_anomaly = gz[:, :, slice_idx_z]
    p = plot(layout = (1, 2), size = (1000, 500), dpi=300)
    heatmap!(p[1, 1], density_map, title="Density Model (Horizontal Slice)", xlabel="X", ylabel="Y", color=:viridis, grid=:true)
    heatmap!(p[1, 2], bouguer_anomaly, title="Bouguer Anomaly (gz)", xlabel="X", ylabel="Y", color=:plasma, grid=:true)
    savefig(p, "PlotsinOne.pdf")
    println("Everything Plotted.pdf.")
end

function GravityForward(filename::String = "gravity_model.ubc")
    model = load_gravity_model_ubc(filename)
    φ = solve_poisson_direct(model)
    gx, gy, gz = compute_gravity(model, φ)
    plot_and_save_results(model, gz)
    return φ, gx, gy, gz
end

φ, gx, gy, gz = GravityForward("gravity_model.ubc")
