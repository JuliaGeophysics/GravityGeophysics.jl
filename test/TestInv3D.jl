# Inversion: control files, objective gradient (4.8), synthetic block (LSQR), file workflow.

# 16 × 16 × 8 mesh of 25 m cells with 2 padding cells per side, a block at 50–125 m depth
function block_problem(; noise = 0.02, seed = 3)
    pad = [100.0, 300.0]
    h = vcat(reverse(pad), fill(25.0, 12), pad)
    m = Model(h, h, vcat(fill(25.0, 6), [50.0, 100.0]); origin = [-sum(h) / 2, -sum(h) / 2, 0.0])
    true_A = zeros(size(m.A)); true_A[7:10, 7:10, 3:5] .= 300.0
    x = range(-140, 140; length = 13)
    d = station_data([a for a in x for _ in x], [b for _ in x for b in x], fill(-20.0, 169))
    op = GravForward(m, d)
    clean = forward(op, vec(true_A))
    σ = noise * maximum(abs, clean) .+ 0.0 .* clean
    rng = Xoshiro(seed)
    d.value .= clean .+ σ .* randn(rng, length(clean)); d.error .= σ
    m, d, true_A
end

@testset "control files" begin
    dir = mktempdir()
    f = FwdCtrl3D(solver = :matrixfree)
    @test ReadFwdCtrl3D(WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), f)) == f
    for c in (InvCtrl3D(algorithm = :lsqr, lambda = 2.0, target_rms = 1.2, max_iter = 7, bounds = (-500.0, 500.0),
                        weighting = :depth, smooth_z = 0.5, lsqr_max_iter = 80),
              InvCtrl3D(algorithm = :lsqr, lambda = 1.0, target_rms = 1.0, max_iter = 50, max_step = 50.0,
                        lsqr_tolerance = 1e-5))
        p = WriteInvCtrl3D(joinpath(dir, "InvCtrl"), c)
        @test ReadInvCtrl3D(p) == c
        @test occursin("Initial damping factor lambda", read(p, String))
        @test occursin("Maximum LSQR iterations", read(p, String))
        @test !occursin("Norm", read(p, String))
    end
    p = joinpath(dir, "InvCtrl")
    write(p, "Algorithm : LSQR\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\n")
    @test_throws ErrorException ReadInvCtrl3D(p)                          # missing required key
    write(p, "Algorithm : LSQR\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\nMaximum number of iterations : 5\nSmoothing weight q : 1\n")
    @test_throws ErrorException ReadInvCtrl3D(p)                          # unknown key
    write(p, "Algorithm : LSQR\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\nMaximum number of iterations : 5\nGN damping : 0\n")
    @test_throws ErrorException ReadInvCtrl3D(p)                          # Gauss-Newton is gone
    write(p, "Algorithm : GN\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\nMaximum number of iterations : 5\n")
    @test_throws ErrorException ReadInvCtrl3D(p)
    write(p, "Algorithm : LSQR\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\nMaximum number of iterations : 5\nNorm : L1\n")
    @test_throws ErrorException ReadInvCtrl3D(p)                          # not an LSQR key
    write(p, "Algorithm : LS\nInitial damping factor lambda : 1\nExit search when rms is less than : 1\nMaximum number of iterations : 5\n")
    @test_throws ErrorException ReadInvCtrl3D(p)                          # LS is now LSQR
    write(p, "Forward solver : sparse\n")
    @test_throws ErrorException ReadFwdCtrl3D(p)
end

@testset "4.8 objective gradient and Hessian" begin
    for weighting in (:sensitivity, :depth, :none)
        m, d, _ = block_problem()
        op = GravForward(m, d)
        P = Inv3DProblem(op, Inv3DOptions(; weighting); start = fill(10.0, length(op.cells)))
        P.beta = 0.3
        x = 10 .+ 20 .* abs.(sin.(1:length(P.var))); v = cos.(1:length(x))
        g = dot(GG.gradient(P, x), v)
        fd(h) = (GG.objective(P, x .+ h .* v) - GG.objective(P, x .- h .* v)) / 2h
        report("objective gradient vs central difference, $weighting", abs(fd(1e-2) - g) / abs(g), 1e-6)
        # the least-squares Hessian is symmetric (exact for linear data)
        u = sin.(2 .* (1:length(x)))
        Hv(v) = GG.data_hessian(P, x, v) .+ P.beta .* GG.reg_hessian(P, v)
        report("Hessian symmetry, $weighting", abs(dot(u, Hv(v)) - dot(v, Hv(u))) / abs(dot(u, Hv(v))), 1e-10)
    end
end

@testset "synthetic block: stops at the target rms" begin
    opt = Inv3DOptions(max_iter = 30, bounds = (0.0, 1000.0), verbose = false)
    for (alg, solver) in ((LSQR3DConfig(), :dense), (LSQR3DConfig(), :matrixfree))
        m, d, true_A = block_problem()
        start = Model(m.dx, m.dy, m.dz; origin = m.origin)
        o = opt
        r = Invert3D(start, d; algorithm = alg, options = o, fwd = FwdCtrl3D(; solver))
        @printf("  %-5s %-10s %3d iterations, rms %.4f, %s\n", inv3d_tag(alg), solver, r.history[end].iteration, r.rms, r.reason)
        @test r.reason == :target_rms && r.converged
        @test (1 - o.target_tolerance) * o.target_rms <= r.rms <= o.target_rms
        @test r isa Inv3DResult{typeof(alg)}
        @test minimum(r.model.A) >= 0                        # positivity bounds
        # peak inside the block (smooth inversions smear it with depth); correlation over the core
        I = argmax(r.model.A)
        @test I[1] in 7:10 && I[2] in 7:10 && I[3] in 3:5
        core = (3:14, 3:14, 1:6)
        @test cor(vec(r.model.A[core...]), vec(true_A[core...])) > 0.5
        @test isapprox(r.predicted, forward(GravForward(m, d), model_vector(GravForward(m, d), r.model)); rtol = 1e-10)
    end
    # the shorthands run the same algorithm
    m, d, _ = block_problem()
    start = Model(m.dx, m.dy, m.dz; origin = m.origin)
    a = LSQR3D(start, d; options = opt)
    @test a.model.A == Invert3D(start, d; options = opt).model.A
end

@testset "LSQR" begin
    # the solver: a dense least-squares problem to round-off, with and without column scaling
    rng = Xoshiro(4)
    A, b = randn(rng, 60, 20), randn(rng, 60)
    for D in (ones(20), 0.1 .+ rand(rng, 20))
        s, k = GG._lsqr(v -> A * v, u -> A' * u, b, D, 100, 1e-14)
        report("LSQR vs QR on a 60 × 20 problem ($k iterations)", norm(s .- A \ b) / norm(A \ b), 1e-10)
    end
    # the stacked system: Aᵀ is the adjoint of A, Aᵀb = −∇Φ and AᵀA the Gauss-Newton Hessian
    m, d, _ = block_problem()
    op = GravForward(m, d)
    P = Inv3DProblem(op, Inv3DOptions(); start = fill(10.0, length(op.cells)))
    P.beta = 0.3
    x = 10 .+ 5 .* sin.(1:length(P.var))
    Av, Atu, b = GG._lsqr3d_system(P, x, trues(length(x)))
    v, u = cos.(1:length(x)), sin.(1:length(b))
    report("LSQR system: adjoint", abs(dot(Av(v), u) - dot(v, Atu(u))) / abs(dot(v, Atu(u))), 1e-10)
    g = GG.gradient(P, x)
    report("LSQR system: Aᵀb = −∇Φ", norm(Atu(b) .+ g) / norm(g), 1e-10)
    Hv = GG.data_hessian(P, x, v) .+ P.beta .* GG.reg_hessian(P, v)
    report("LSQR system: AᵀA = GN Hessian", norm(Atu(Av(v)) .- Hv) / norm(Hv), 1e-10)
end

@testset "first step overshooting the target is redone" begin
    m, d, _ = block_problem()
    d.error .*= 4                                            # errors 4× too large: the first LS step lands far below rms 1
    start = Model(m.dx, m.dy, m.dz; origin = m.origin)
    o = Inv3DOptions(verbose = false)
    r = Invert3D(start, d; options = o)
    @test r.reason == :target_rms && (1 - o.target_tolerance) * o.target_rms <= r.rms <= o.target_rms
end

@testset "mask keeps fixed cells" begin
    m, d, _ = block_problem()
    start = Model(m.dx, m.dy, m.dz, fill(20.0, size(m.A)); origin = m.origin)
    mask = trues(size(m.A)); mask[:, :, 1] .= false
    r = Invert3D(start, d; options = Inv3DOptions(max_iter = 5, verbose = false), mask)
    @test all(r.model.A[:, :, 1] .== 20.0)
end

@testset "per-cell bounds (drillhole constraint)" begin
    m, d, _ = block_problem()
    start = Model(m.dx, m.dy, m.dz; origin = m.origin)
    lo, hi = fill(0.0, size(m.A)), fill(1000.0, size(m.A))
    hi[7:10, 7:10, 4:5] .= 10.0                               # the block's lower part held light
    lo[2:3, 2:3, 1:2] .= 200.0                                # off the block, held dense
    lower, upper = (Model(m.dx, m.dy, m.dz, A; origin = m.origin) for A in (lo, hi))
    r = Invert3D(start, d; options = Inv3DOptions(max_iter = 10, verbose = false), lower, upper)
    @test all(r.model.A[7:10, 7:10, 4:5] .<= 10.0) && all(r.model.A[2:3, 2:3, 1:2] .>= 200.0)
    @test all(lo .<= r.model.A .<= hi)
    @test_throws ArgumentError Invert3D(start, d; options = Inv3DOptions(verbose = false), lower = upper, upper = lower)
end

@testset "file workflow" begin
    m, d, _ = block_problem()
    dir = mktempdir()
    sp = write_model(joinpath(dir, "model.start"), Model(m.dx, m.dy, m.dz; origin = m.origin))
    pp = write_model(joinpath(dir, "model.prior"), Model(m.dx, m.dy, m.dz; origin = m.origin))
    mp = write_mask(joinpath(dir, "model.mask"), m, active_cells(m))
    dp = write_data_grav(joinpath(dir, "obs.dat"), d)
    fp = WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), FwdCtrl3D())
    cp = WriteInvCtrl3D(joinpath(dir, "InvCtrl.LSQR"), InvCtrl3D(algorithm = :lsqr, lambda = 1.0, target_rms = 1.0, max_iter = 20))
    r = Invert3D(sp, dp, fp, cp, mp, pp; run_dir = joinpath(dir, "run"))
    for f in ("model.den", "data.pred", "0lsqr.log", "inputs/model.start", "inputs/obs.dat", "inputs/FwdCtrl",
              "inputs/InvCtrl.LSQR", "inputs/model.mask", "inputs/model.prior")
        @test isfile(joinpath(r.run_dir, f))
    end
    mod = load_model(joinpath(r.run_dir, "model.den"))
    @test same_mesh(mod, m) && maximum(abs, mod.A .- r.final.A) < 1e-9 * maximum(abs, r.final.A)
    pred = load_data_grav(joinpath(r.run_dir, "data.pred"))
    @test isapprox(pred.value, r.predicted.value; rtol = 1e-5)
    @test r.algorithm == :lsqr && r.converged && r.reason == :target_rms
    @test !isfile(joinpath(r.run_dir, "History.csv")) && !isfile(joinpath(r.run_dir, "Summary.txt"))
    lines = readlines(joinpath(r.run_dir, "0lsqr.log"))
    rows = filter(l -> tryparse(Int, first(split(l))) !== nothing, lines)
    @test startswith(lines[1], "# LSQR 3D gravity inversion") && length(rows) == length(r.history)
    @test parse(Float64, split(rows[end])[4]) ≈ r.rms rtol = 1e-5
    @test any(l -> occursin(r"^# termination\s+= target_rms$", l), lines) && any(startswith("# rms "), lines)
    # optional bounds models, copied with the inputs
    hi = fill(1000.0, size(m.A)); hi[:, :, 1] .= 5.0
    lp = write_model(joinpath(dir, "model.lower"), Model(m.dx, m.dy, m.dz; origin = m.origin))
    up = write_model(joinpath(dir, "model.upper"), Model(m.dx, m.dy, m.dz, hi; origin = m.origin))
    rb = Invert3D(sp, dp, fp, cp, mp, pp; run_dir = joinpath(dir, "runb"), bounds = (lp, up))
    @test isfile(joinpath(rb.run_dir, "inputs/model.lower")) && isfile(joinpath(rb.run_dir, "inputs/model.upper"))
    @test all(rb.final.A[:, :, 1] .<= 5.0)
    bad = write_model(joinpath(dir, "other.prior"), Model(m.dx[2:end], m.dy, m.dz; origin = m.origin))
    @test_throws ErrorException Invert3D(sp, dp, fp, cp, mp, bad; run_dir = joinpath(dir, "run2"))
end
