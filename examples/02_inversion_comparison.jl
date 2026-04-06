using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using GravityGeophysics
using Printf

function main()
    outdir = joinpath(@__DIR__, "demo_output", "dipping_block_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    mkpath(outdir)

    bundle = generate_block_benchmark_bundle(outdir; noise_fraction = 0.01, seed = 42)
    true_model = bundle[:model]
    observed = bundle[:data]

    G = sensitivity_matrix_prism(true_model, observed.receivers)
    standard_model, inr_model, comparison = compare_inversion_methods(
        observed,
        true_model;
        G = G,
        tikhonov_kwargs = (lambda = 0.25, alpha_s = 0.02, alpha_xyz = (1.0, 1.0, 1.5), depth_beta = 1.0),
        inr_kwargs = (epochs = 500, learning_rate = 1e-2, hidden_dim = 256, depth = 4, nfreq = 2, max_density = 600.0, seed = 42, verbose = false),
    )

    save_model_vox(joinpath(outdir, "true_density_contrast.vox"), true_model; name = "True Dipping Block Density Contrast", kind = :density_contrast)
    save_model_vox(joinpath(outdir, "standard_inversion_contrast.vox"), standard_model; name = "Standard Tikhonov Density Contrast", kind = :density_contrast)
    save_model_vox(joinpath(outdir, "inr_inversion_contrast.vox"), inr_model; name = "INR Density Contrast", kind = :density_contrast)

    z_slice_index = min(true_model.mesh.nz, 6)
    x_slice_index = div(true_model.mesh.nx, 2) + 1
    y_slice_index = div(true_model.mesh.ny, 2) + 1

    plot_orthogonal_inversion_comparison(
        true_model,
        standard_model,
        inr_model;
        x_slice_index = x_slice_index,
        y_slice_index = y_slice_index,
        z_slice_index = z_slice_index,
        save_path = joinpath(outdir, "orthogonal_inversion_comparison.png"),
        use_bulk_density = false,
        clims = (0.0, 400.0),
    )
    plot_inversion_comparison(true_model, standard_model, inr_model; slice_axis = :z, slice_index = z_slice_index, save_path = joinpath(outdir, "inversion_comparison.png"))
    plot_convergence(comparison["inr"]; save_path = joinpath(outdir, "inr_convergence.png"))
    plot_gravity_data(observed; title = "Observed Gravity, Dipping Block Benchmark", save_path = joinpath(outdir, "observed_gravity.png"))

    open(joinpath(outdir, "inversion_summary.txt"), "w") do io
        println(io, "Output directory: " * outdir)
        println(io, "Scenario: Python-style dipping staircase block benchmark")
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