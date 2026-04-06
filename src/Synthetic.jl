function create_terrascope_demo_mesh(; nx::Int = 24, ny::Int = 24, nz::Int = 16,
    lx::Float64 = 160000.0, ly::Float64 = 160000.0, lz::Float64 = 50000.0,
    x0::Float64 = -80000.0, y0::Float64 = -80000.0, z0::Float64 = 0.0)
    return create_mesh(nx = nx, ny = ny, nz = nz, lx = lx, ly = ly, lz = lz, x0 = x0, y0 = y0, z0 = z0)
end

function create_block_benchmark_mesh(; dx::Float64 = 50.0, dy::Float64 = 50.0, dz::Float64 = 50.0,
    lx::Float64 = 1000.0, ly::Float64 = 1000.0, lz::Float64 = 500.0,
    x_center_min::Float64 = 385000.0, y_center_min::Float64 = 6670000.0, z_center_min::Float64 = 0.0)
    nx = Int(round(lx / dx)) + 1
    ny = Int(round(ly / dy)) + 1
    nz = Int(round(lz / dz)) + 1
    dxv = fill(dx, nx)
    dyv = fill(dy, ny)
    dzv = fill(dz, nz)
    return create_mesh(
        nx = nx,
        ny = ny,
        nz = nz,
        dx = dxv,
        dy = dyv,
        dz = dzv,
        x0 = x_center_min - 0.5 * dx,
        y0 = y_center_min - 0.5 * dy,
        z0 = z_center_min - 0.5 * dz,
    )
end

function create_block_benchmark_model(mesh::GravityMesh = create_block_benchmark_mesh(); background::Float64 = 0.0, block_density::Float64 = 400.0)
    model = create_model(mesh; background = background)
    for step in 0:6
        k = 2 + step
        y_start = 12 - step
        y_end = 16 - step
        x_start = 8
        x_end = 13
        if 1 <= k <= mesh.nz
            for j in max(1, y_start):min(mesh.ny, y_end), i in max(1, x_start):min(mesh.nx, x_end)
                model.density[i, j, k] = block_density
            end
        end
    end
    return model
end

function create_block_benchmark_receivers(mesh::GravityMesh; z::Float64 = -1.0)
    xc, yc, _ = get_cell_centers(mesh)
    x = repeat(xc, inner = length(yc))
    y = repeat(yc, outer = length(xc))
    zvals = fill(z, length(x))
    return create_receivers(x, y, zvals)
end

function create_pyhasalmi_mesh(; dx::Float64 = 100.0, dy::Float64 = 100.0, dz::Float64 = 80.0,
    lx::Float64 = 2800.0, ly::Float64 = 2800.0, lz::Float64 = 1600.0,
    x_center::Float64 = 452528.6291993095, y_center::Float64 = 7059329.16324865, z0::Float64 = 0.0)
    nx = Int(round(lx / dx))
    ny = Int(round(ly / dy))
    nz = Int(round(lz / dz))
    return create_mesh(
        nx = nx,
        ny = ny,
        nz = nz,
        dx = fill(dx, nx),
        dy = fill(dy, ny),
        dz = fill(dz, nz),
        x0 = x_center - 0.5 * lx,
        y0 = y_center - 0.5 * ly,
        z0 = z0,
    )
end

function create_pyhasalmi_synthetic_model(mesh::GravityMesh = create_pyhasalmi_mesh(); background::Float64 = 0.0)
    model = create_model(mesh; background = background)
    xc, yc, zc = get_cell_centers(mesh)
    x_center = 452528.6291993095
    y_center = 7059329.16324865
    strike = deg2rad(32.0)
    cos_strike = cos(strike)
    sin_strike = sin(strike)

    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        x_local = xc[i] - x_center
        y_local = yc[j] - y_center
        z_local = zc[k]

        along = cos_strike * x_local + sin_strike * y_local
        cross = -sin_strike * x_local + cos_strike * y_local

        ore_lens = 620.0 * exp(-0.5 * ((along - 120.0) / 290.0)^2 - 0.5 * (cross / 110.0)^2 - 0.5 * ((z_local - 620.0) / 170.0)^2)
        deep_pod = 360.0 * exp(-0.5 * ((along + 180.0) / 180.0)^2 - 0.5 * ((cross - 20.0) / 90.0)^2 - 0.5 * ((z_local - 980.0) / 180.0)^2)
        alteration_halo = 120.0 * exp(-0.5 * (along / 620.0)^2 - 0.5 * (cross / 260.0)^2 - 0.5 * ((z_local - 720.0) / 320.0)^2)
        hangingwall_low = -95.0 * exp(-0.5 * ((along + 40.0) / 420.0)^2 - 0.5 * ((cross + 180.0) / 150.0)^2 - 0.5 * ((z_local - 260.0) / 120.0)^2)
        shallow_overburden = -35.0 * exp(-0.5 * (x_local / 900.0)^2 - 0.5 * (y_local / 900.0)^2 - 0.5 * ((z_local - 60.0) / 50.0)^2)

        model.density[i, j, k] += ore_lens + deep_pod + alteration_halo + hangingwall_low + shallow_overburden
    end

    return model
end

function create_pyhasalmi_receivers(mesh::GravityMesh; nrx::Int = 31, nry::Int = 31, z::Float64 = -5.0, xfrac::Float64 = 0.82, yfrac::Float64 = 0.82)
    return create_receivers(mesh; nrx = nrx, nry = nry, z = z, xfrac = xfrac, yfrac = yfrac)
end

function create_terrascope_synthetic_model(mesh::GravityMesh = create_terrascope_demo_mesh(); background::Float64 = 0.0)
    model = create_model(mesh; background = background)

    add_block_anomaly!(model;
        xmin = -18000.0, xmax = 12000.0,
        ymin = -26000.0, ymax = -2000.0,
        zmin = 4000.0, zmax = 18000.0,
        drho = 320.0)

    add_block_anomaly!(model;
        xmin = 18000.0, xmax = 42000.0,
        ymin = 10000.0, ymax = 30000.0,
        zmin = 12000.0, zmax = 28000.0,
        drho = -210.0)

    add_sphere_anomaly!(model;
        cx = -35000.0, cy = 30000.0, cz = 26000.0,
        radius = 9000.0,
        drho = 260.0)

    add_sphere_anomaly!(model;
        cx = 8000.0, cy = -12000.0, cz = 34000.0,
        radius = 11000.0,
        drho = 180.0)

    xc, yc, zc = get_cell_centers(mesh)
    for k in 1:mesh.nz, j in 1:mesh.ny, i in 1:mesh.nx
        x = xc[i]
        y = yc[j]
        z = zc[k]
        model.density[i, j, k] += 45.0 * exp(-((x - 26000.0)^2) / (2 * (18000.0^2))) * exp(-((y + 5000.0)^2) / (2 * (12000.0^2))) * (z / maximum(zc))
    end

    return model
end

function create_demo_receivers(mesh::GravityMesh; nrx::Int = 41, nry::Int = 41, z::Float64 = -50.0, xfrac::Float64 = 0.92, yfrac::Float64 = 0.92)
    return create_receivers(mesh; nrx = nrx, nry = nry, z = z, xfrac = xfrac, yfrac = yfrac)
end

function generate_demo_bundle(outdir::AbstractString; noise_std_mgal::Float64 = 0.03)
    mkpath(outdir)
    mesh = create_terrascope_demo_mesh()
    model = create_terrascope_synthetic_model(mesh)
    receivers = create_demo_receivers(mesh)
    clean = forward(model, receivers, Prism())
    noise_std_si = noise_std_mgal / SI_TO_MGAL
    noisy = GravityData(receivers, clean.gz .+ noise_std_si .* randn(length(clean.gz)), fill(noise_std_si, length(clean.gz)))

    model_vox_path = joinpath(outdir, "terrascope_demo_model.vox")
    data_xyz_path = joinpath(outdir, "terrascope_demo_data.xyz")
    recv_xyz_path = joinpath(outdir, "terrascope_demo_receivers.xyz")

    save_model_vox(model_vox_path, model; name = "TerraScope Gravity Synthetic", kind = :density_contrast)
    save_data_xyz(data_xyz_path, noisy)
    save_receivers_xyz(recv_xyz_path, receivers)

    return Dict(
        :mesh => mesh,
        :model => model,
        :receivers => receivers,
        :data => noisy,
        :model_vox_path => model_vox_path,
        :data_xyz_path => data_xyz_path,
        :receivers_xyz_path => recv_xyz_path,
    )
end

function generate_block_benchmark_bundle(outdir::AbstractString; noise_fraction::Float64 = 0.01, seed::Int = 42)
    mkpath(outdir)
    Random.seed!(seed)

    mesh = create_block_benchmark_mesh()
    model = create_block_benchmark_model(mesh)
    receivers = create_block_benchmark_receivers(mesh)
    clean = forward(model, receivers, Prism())
    noise_std_si = noise_fraction * std(clean.gz)
    noisy = GravityData(receivers, clean.gz .+ noise_std_si .* randn(length(clean.gz)), fill(noise_std_si, length(clean.gz)))

    model_vox_path = joinpath(outdir, "block_benchmark_model.vox")
    data_xyz_path = joinpath(outdir, "block_benchmark_data.xyz")
    recv_xyz_path = joinpath(outdir, "block_benchmark_receivers.xyz")

    save_model_vox(model_vox_path, model; name = "INR Gravity Benchmark", kind = :density_contrast)
    save_data_xyz(data_xyz_path, noisy)
    save_receivers_xyz(recv_xyz_path, receivers)

    return Dict(
        :mesh => mesh,
        :model => model,
        :receivers => receivers,
        :data => noisy,
        :model_vox_path => model_vox_path,
        :data_xyz_path => data_xyz_path,
        :receivers_xyz_path => recv_xyz_path,
        :noise_std_si => noise_std_si,
        :noise_fraction => noise_fraction,
    )
end

function generate_pyhasalmi_bundle(outdir::AbstractString; noise_std_mgal::Float64 = 0.025, seed::Int = 42,
    mesh::GravityMesh = create_pyhasalmi_mesh(), receiver_kwargs::NamedTuple = (;), host_density::Float64 = 2670.0)
    mkpath(outdir)
    Random.seed!(seed)

    model = create_pyhasalmi_synthetic_model(mesh)
    receivers = create_pyhasalmi_receivers(mesh; pairs(receiver_kwargs)...)
    clean = forward(model, receivers, Prism())
    noise_std_si = noise_std_mgal / SI_TO_MGAL
    noisy = GravityData(receivers, clean.gz .+ noise_std_si .* randn(length(clean.gz)), fill(noise_std_si, length(clean.gz)))

    model_vox_path = joinpath(outdir, "pyhasalmi_true_density_contrast.vox")
    data_xyz_path = joinpath(outdir, "pyhasalmi_observed_gravity.xyz")
    recv_xyz_path = joinpath(outdir, "pyhasalmi_receivers.xyz")

    save_model_vox(model_vox_path, model; name = "Pyhasalmi synthetic density contrast", kind = :density_contrast)
    save_data_xyz(data_xyz_path, noisy)
    save_receivers_xyz(recv_xyz_path, receivers)

    return Dict(
        :mesh => mesh,
        :model => model,
        :receivers => receivers,
        :data => noisy,
        :model_vox_path => model_vox_path,
        :data_xyz_path => data_xyz_path,
        :receivers_xyz_path => recv_xyz_path,
        :noise_std_si => noise_std_si,
        :noise_std_mgal => noise_std_mgal,
        :host_density => host_density,
        :center_easting => 452528.6291993095,
        :center_northing => 7059329.16324865,
    )
end