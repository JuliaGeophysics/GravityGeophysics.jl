# Control files
# Author: @pankajkmishra
# Round trips of FwdCtrl, InvCtrl.GN, InvCtrl.INR and the covariance, the ReadInvCtrl dispatch on the
# Algorithm key, the defaults of absent keys, and the validation of bad values

using Test

@testset "control files" begin
    mktempdir() do dir
        fwd = FwdCtrl3D(write_frechet = true)
        @test ReadFwdCtrl3D(WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), fwd)) == fwd

        gn = InvCtrl3D(algorithm = :gn, lambda = 12.5, target_rms = 1.2, max_iter = 7, lambda_cooling = 3.0,
                       bounds = (-250.0, 750.0), smallness = 0.04, smooth_x = 2.0, smooth_y = 0.5, smooth_z = 1.5,
                       weighting = :depth, depth_exponent = 3.0, depth_z0 = 25.0, error_floor = 0.02,
                       error_percent = 3.0, gn_damping = 1e-2, cg_iterations = 30, cg_tolerance = 1e-4)
        @test ReadInvCtrl3D(WriteInvCtrl3D(joinpath(dir, "InvCtrl.GN"), gn)) == gn

        inr = InvCtrlINR3D(algorithm = :inr, target_rms = 0.9, max_iter = 300, bounds = (-400.0, 400.0),
                           error_floor = 0.05, error_percent = 1.0, width = 64, layers = 3, frequencies = 4,
                           learning_rate = 2e-3, scaling = :isotropic, seed = 7, divergence_factor = 2.5,
                           max_restarts = 1, schedule = :none, rate_floor = 0.1, anneal = 0.0)
        @test ReadInvCtrlINR3D(WriteInvCtrlINR3D(joinpath(dir, "InvCtrl.INR"), inr)) == inr

        # the training keys default when an older file leaves them out
        older = replace(read(joinpath(dir, "InvCtrl.INR"), String),
                        r"(?m)^(Learning rate schedule|Learning rate floor|Frequency annealing).*\n" => "")
        write(joinpath(dir, "old.INR"), older)
        plain = ReadInvCtrlINR3D(joinpath(dir, "old.INR"))
        @test (plain.schedule, plain.rate_floor, plain.anneal) == (:cosine, 0.05, 0.5)

        # the dispatcher picks the type from the Algorithm key
        @test ReadInvCtrl(joinpath(dir, "InvCtrl.GN")) isa InvCtrl3D
        @test ReadInvCtrl(joinpath(dir, "InvCtrl.INR")) isa InvCtrlINR3D

        # lambda auto, and keys left out take their defaults
        auto = InvCtrl3D(algorithm = :gn, lambda = nothing, target_rms = 1.0, max_iter = 5)
        back = ReadInvCtrl3D(WriteInvCtrl3D(joinpath(dir, "auto.GN"), auto))
        @test back.lambda === nothing && back.weighting == :sensitivity && back.lambda_cooling == 2.0
        open(joinpath(dir, "min.GN"), "w") do io
            println(io, "Algorithm : GN")
            println(io, "Initial damping factor lambda : auto")
            println(io, "Exit search when rms is less than : 1.0")
            println(io, "Maximum number of iterations : 5")
        end
        @test ReadInvCtrl3D(joinpath(dir, "min.GN")) == auto

        # 'Invert padding' round trips on both algorithms and defaults true for GN, false for the INR
        @test auto.invert_padding && !ReadInvCtrlINR3D(joinpath(dir, "InvCtrl.INR")).invert_padding
        @test ReadInvCtrl3D(joinpath(dir, "min.GN")).invert_padding
        gn_fixed = InvCtrl3D(algorithm = :gn, lambda = nothing, target_rms = 1.0, max_iter = 5,
                             invert_padding = false)
        @test !ReadInvCtrl3D(WriteInvCtrl3D(joinpath(dir, "pad.GN"), gn_fixed)).invert_padding
        inr_free = InvCtrlINR3D(algorithm = :inr, target_rms = 1.0, max_iter = 5, invert_padding = true)
        @test ReadInvCtrlINR3D(WriteInvCtrlINR3D(joinpath(dir, "pad.INR"), inr_free)).invert_padding

        # a control file written before the key existed reads with the per-algorithm default
        for (src, dst, reader, want) in (("InvCtrl.GN", "nopad.GN", ReadInvCtrl3D, true),
                                         ("InvCtrl.INR", "nopad.INR", ReadInvCtrlINR3D, false))
            write(joinpath(dir, dst), replace(read(joinpath(dir, src), String), r"(?m)^Invert padding.*\n" => ""))
            @test reader(joinpath(dir, dst)).invert_padding == want
        end

        # covariance round trip, including a layer-varying mask
        mask = ones(Int, 4, 5, 3)
        mask[:, :, 1] .= 0
        mask[2, 3, 2] = 2
        cov = Cov3D(mask; smoothing = 0.25, n_smooth = 3)
        back_cov = ReadCov3D(WriteCov3D(joinpath(dir, "C3.dat"), cov))
        @test back_cov.mask == mask && back_cov.sz == 0.25 && back_cov.n_smooth == 3

        # validation
        @test_throws ArgumentError InvCtrl3D(algorithm = :inr, lambda = nothing, target_rms = 1.0, max_iter = 2) |>
                                   c -> WriteInvCtrl3D(joinpath(dir, "bad.GN"), c)
        @test_throws ArgumentError WriteInvCtrl3D(joinpath(dir, "bad2.GN"),
                                                  InvCtrl3D(algorithm = :gn, lambda = nothing, target_rms = 1.0,
                                                            max_iter = 2, bounds = (10.0, -10.0)))
        @test_throws ArgumentError WriteInvCtrlINR3D(joinpath(dir, "bad.INR"),
                                                     InvCtrlINR3D(algorithm = :inr, target_rms = 1.0, max_iter = 2,
                                                                  divergence_factor = 0.5))
        @test_throws ErrorException ReadFwdCtrl3D(joinpath(dir, "missing"))
    end
end
