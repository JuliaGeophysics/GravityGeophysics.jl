using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using GravityGeophysics
using Printf

function main()
    outdir = joinpath(@__DIR__, "demo_output", "forward_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    mkpath(outdir)

    bundle = generate_demo_bundle(outdir)
    mesh = bundle[:mesh]
    model = bundle[:model]
    data_obs = bundle[:data]
    receivers = bundle[:receivers]

    data_prism = forward(model, receivers, Prism())
    data_fd, fd_info = forward_fd(model, receivers; options = FDOptions(tol = 1e-10, maxiter = 6000, use_prism_bc = false))
    stats = compute_stats(data_prism, data_fd)

    prism_obs = GravityData(receivers, data_prism.gz, fill(0.0, receivers.n))
    fd_obs = GravityData(receivers, data_fd.gz, fill(0.0, receivers.n))

    save_data_xyz(joinpath(outdir, "gravity_prism.xyz"), prism_obs)
    save_data_xyz(joinpath(outdir, "gravity_fd.xyz"), fd_obs)

    plot_gravity_comparison(prism_obs, fd_obs; title_prefix = "TerraScope-area Gravity", save_path = joinpath(outdir, "gravity_comparison.png"))
    plot_model_slice(model; slice_axis = :z, slice_index = div(mesh.nz, 3), save_path = joinpath(outdir, "model_slice_z.png"))
    plot_gravity_data(data_obs; title = "Synthetic Observed Gravity", save_path = joinpath(outdir, "synthetic_observed_gravity.png"))
    plot_receiver_locations(receivers; save_path = joinpath(outdir, "receiver_locations.png"))

    open(joinpath(outdir, "forward_summary.txt"), "w") do io
        println(io, "Output directory: " * outdir)
        @printf(io, "Prism max response: %.6f mGal\n", maximum(data_prism.gz) * SI_TO_MGAL)
        @printf(io, "FD max response: %.6f mGal\n", maximum(data_fd.gz) * SI_TO_MGAL)
        @printf(io, "FD iterations: %d\n", fd_info["iterations"])
        @printf(io, "RMSE (Prism vs FD): %.6f mGal\n", stats["rmse_mgal"])
        @printf(io, "Relative RMSE: %.4f %%\n", stats["rel_rmse_pct"])
    end

    println("Saved forward example outputs to: " * outdir)
    println("RMSE Prism vs FD: $(round(stats["rmse_mgal"], digits = 6)) mGal")
    return outdir, stats, fd_info
end

main()
