# 3D file workflows
# Author: @pankajkmishra
# Runs ForwardSolve3D and the six-file Invert3D on a small case, once per algorithm: the run directory holds the
# inputs, model.den, data.pred, History.log and Summary.txt, and the covariance mask fixes cells

using Test

@testset "3D file workflows" begin
    mktempdir() do dir
        truth = small_mesh(n = 8, nz = 5)
        truth.A[3:5, 3:5, 2:3] .= 300.0
        d = small_survey(n = 9, h = 45.0)
        d.value .= forward3d(truth, d)
        d.error .= 0.02
        start = WriteModel3D(joinpath(dir, "model.start"), model_like(truth, zeros(size(truth.A))))
        data = WriteData3D(joinpath(dir, "data.xyz"), d)
        fwd = WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), FwdCtrl3D(write_frechet = true))
        inv = WriteInvCtrl3D(joinpath(dir, "InvCtrl.GN"),
                             InvCtrl3D(algorithm = :gn, lambda = nothing, target_rms = 1.0, max_iter = 15,
                                       bounds = (-500.0, 500.0)))
        true_path = WriteModel3D(joinpath(dir, "model.true"), truth)
        mask = ones(Int, size(truth.A))
        mask[:, :, end] .= 0
        cov = WriteCov3D(joinpath(dir, "C3.dat"), Cov3D(mask))

        pred = ForwardSolve3D(true_path, data, fwd)
        @test ReadData3D(pred).value ≈ d.value rtol = 1e-4
        @test LoadSensitivity3D(joinpath(dir, "sensitivity.jld2")).sensitivity |> size == (length(truth.A), ndata(d))

        run = Invert3D(start, data, fwd, inv, cov, start; run_dir = joinpath(dir, "run"))
        @test run.algorithm == :gn && run.rms <= 1.0 && run.converged
        for f in ("model.den", "data.pred", "History.log", "Summary.txt", "sensitivity.jld2",
                  joinpath("inputs", "FwdCtrl"), joinpath("inputs", "C3.dat"))
            @test isfile(joinpath(run.run_dir, f))
        end
        @test all(ReadModel3D(joinpath(run.run_dir, "model.den")).A[:, :, end] .== 0.0)
        @test occursin("Termination: target_rms", read(joinpath(run.run_dir, "Summary.txt"), String))

        # the same six files with an InvCtrl.INR are dispatched to the INR driver
        inr = WriteInvCtrlINR3D(joinpath(dir, "InvCtrl.INR"),
                                InvCtrlINR3D(algorithm = :inr, target_rms = 1.0, max_iter = 300,
                                             bounds = (-500.0, 500.0), width = 64, layers = 3,
                                             learning_rate = 1e-2))
        run_inr = Invert3D(start, data, fwd, inr, cov, start; run_dir = joinpath(dir, "run_inr"))
        @test run_inr.algorithm == :inr && run_inr.rms < run_inr.history[1].rms
        @test isfile(joinpath(run_inr.run_dir, "model.den")) && isfile(joinpath(run_inr.run_dir, "data.pred"))
        @test occursin("Algorithm: INR", read(joinpath(run_inr.run_dir, "Summary.txt"), String))
        @test all(ReadModel3D(joinpath(run_inr.run_dir, "model.den")).A[:, :, end] .== 0.0)

        # the default run directory sits next to the data
        @test startswith(basename(GravityGeophysics._inv3d_run_dir(data)), "run_")
    end
end
