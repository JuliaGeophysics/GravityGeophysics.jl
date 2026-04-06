using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using GravityGeophysics
using Printf

function add_density_offset(model::GravityModel, rho0::Float64)
    return GravityModel(model.mesh, model.density .+ rho0)
end

function strongest_depth_slice(model::GravityModel)
    amplitudes = [maximum(abs.(model.density[:, :, k])) for k in 1:model.mesh.nz]
    return argmax(amplitudes)
end

function main()
    outdir = joinpath(@__DIR__, "demo_output", "pyhasalmi_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    mkpath(outdir)

    bundle = generate_pyhasalmi_bundle(outdir; noise_std_mgal = 0.025, seed = 42)
    true_model = bundle[:model]
    observed = bundle[:data]
    host_density = bundle[:host_density]

    G = sensitivity_matrix_prism(true_model, observed.receivers)
    standard_model, inr_model, comparison = compare_inversion_methods(
        observed,
        true_model;
        G = G,
        tikhonov_kwargs = (lambda = 0.15, alpha_s = 0.01, alpha_xyz = (1.0, 1.0, 1.25), depth_beta = 1.0),
        inr_kwargs = (epochs = 650, learning_rate = 5e-3, hidden_dim = 256, depth = 4, nfreq = 2, max_density = 750.0, seed = 42, verbose = false),
    )

    true_bulk = add_density_offset(true_model, host_density)
    standard_bulk = add_density_offset(standard_model, host_density)
    inr_bulk = add_density_offset(inr_model, host_density)

    save_model_vox(joinpath(outdir, "true_density_contrast.vox"), true_model; name = "True Synthetic Density Contrast", kind = :density_contrast)
    save_model_vox(joinpath(outdir, "standard_inversion_contrast.vox"), standard_model; name = "Standard Tikhonov Density Contrast", kind = :density_contrast)
    save_model_vox(joinpath(outdir, "inr_inversion_contrast.vox"), inr_model; name = "INR Density Contrast", kind = :density_contrast)
    save_model_vox(joinpath(outdir, "true_bulk_density.vox"), true_bulk; name = "True Bulk Density", kind = :density)
    save_model_vox(joinpath(outdir, "standard_inversion_bulk_density.vox"), standard_bulk; name = "Standard Tikhonov Bulk Density", kind = :density)
    save_model_vox(joinpath(outdir, "inr_inversion_bulk_density.vox"), inr_bulk; name = "INR Bulk Density", kind = :density)

    z_slice_index = strongest_depth_slice(true_model)
    x_slice_index = div(true_model.mesh.nx, 2) + 1
    y_slice_index = div(true_model.mesh.ny, 2) + 1

    plot_orthogonal_inversion_comparison(
        true_bulk,
        standard_bulk,
        inr_bulk;
        x_slice_index = x_slice_index,
        y_slice_index = y_slice_index,
        z_slice_index = z_slice_index,
        save_path = joinpath(outdir, "orthogonal_inversion_comparison.png"),
        use_bulk_density = true,
    )
    plot_inversion_comparison(true_bulk, standard_bulk, inr_bulk; slice_axis = :z, slice_index = z_slice_index, save_path = joinpath(outdir, "inversion_comparison.png"), use_bulk_density = true)
    plot_convergence(comparison["inr"]; save_path = joinpath(outdir, "inr_convergence.png"))
    plot_gravity_data(observed; title = "Observed Gravity, Pyhasalmi Synthetic", save_path = joinpath(outdir, "observed_gravity.png"))

    open(joinpath(outdir, "inversion_summary.txt"), "w") do io
        println(io, "Output directory: " * outdir)
        println(io, "Scenario: Pyhasalmi-area synthetic sulfide-style anomaly in TM35FIN")
        @printf(io, "Reference center: E=%.3f m, N=%.3f m\n", bundle[:center_easting], bundle[:center_northing])
        @printf(io, "Display host density: %.1f kg/m^3\n", host_density)
        @printf(io, "Slices: z=%d, y=%d, x=%d\n", z_slice_index, y_slice_index, x_slice_index)
        @printf(io, "Standard model RMSE: %.4f kg/m^3\n", comparison["standard_rmse_model"])
        @printf(io, "INR model RMSE: %.4f kg/m^3\n", comparison["inr_rmse_model"])
        @printf(io, "Standard data RMSE: %.6f mGal\n", comparison["standard_rmse_data_mgal"])
        @printf(io, "INR data RMSE: %.6f mGal\n", comparison["inr_rmse_data_mgal"])
    end

    println("Saved inversion outputs to: " * outdir)
    println("Standard model RMSE: $(round(comparison["standard_rmse_model"], digits = 4)) kg/m^3")
    println("INR model RMSE: $(round(comparison["inr_rmse_model"], digits = 4)) kg/m^3")
    return outdir, comparison
end

main()