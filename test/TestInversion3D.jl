# 3D GN and INR inversion
# Author: @pankajkmishra
# Inverts noise-free data of a buried block on a small mesh: the target misfit is reached, the bounds hold, the
# anomaly lands in the right place, fixed cells keep their values, and the gradient matches finite differences
# The INR set checks the same block with the coordinate network: the misfit drops and the body is laterally right
# (depth is not asserted -- without depth weighting the INR pulls mass shallow, which is expected)
# A second INR set inverts noisy data of two blocks of opposite sign: the network stops at the discrepancy
# principle, recovers both signs in the right corners and matches Gauss-Newton's correlation with the truth

using Test
using LinearAlgebra
using Random
using SparseArrays
using Statistics

@testset "3D GN inversion" begin
    truth, d = block_case()
    start = model_like(truth, zeros(size(truth.A)))

    opt = Inv3DOptions(max_iter = 25, target_rms = 1.0, bounds = (-500.0, 500.0), verbose = false)
    r = Invert3D(start, d, PRISM; options = opt)
    @test r.converged && r.reason == :target_rms && r.rms <= 1.0
    @test all(-500 .<= r.density .<= 500)
    @test r.history[1].rms > 3 && issorted([h.beta for h in r.history]; rev = true)
    peak = Tuple(CartesianIndices(r.density)[argmax(r.density)])
    @test 3 <= peak[1] <= 7 && 4 <= peak[2] <= 8
    @test cor(vec(r.density), vec(truth.A)) > 0.5

    # the reported rms is the misfit of the reported predicted data
    @test chi2_rms3d(d, data_like(d, r.predicted)).rms ≈ r.rms rtol = 1e-6

    # bounds tighter than the truth are honoured
    rb = Invert3D(start, d, PRISM; options = Inv3DOptions(max_iter = 6, bounds = (0.0, 100.0), verbose = false))
    @test maximum(rb.density) <= 100 + 1e-12 && minimum(rb.density) >= -1e-12

    # depth weighting runs and a fixed top layer keeps its value
    free = trues(size(truth.A))
    free[:, :, 1] .= false
    s2 = model_like(start, ifelse.(free, 0.0, 25.0))
    rd = Invert3D(s2, d, PRISM; active_cells = free,
                  options = Inv3DOptions(max_iter = 25, target_rms = 1.0, bounds = (-500.0, 500.0),
                                         weighting = :depth, verbose = false))
    @test all(rd.density[:, :, 1] .== 25.0) && rd.rms <= 1.0

    # air cells are excluded and the data of the fixed earth cells is taken off once
    topo = model_like(truth, copy(truth.A))
    topo.A[1:3, 1:3, 1] .= NaN
    rt = Invert3D(model_like(topo, ifelse.(isnan.(topo.A), NaN, 0.0)), d, PRISM;
                  options = Inv3DOptions(max_iter = 15, target_rms = 1.0, bounds = (-500.0, 500.0), verbose = false))
    @test all(isnan, rt.density[1:3, 1:3, 1]) && rt.rms <= 1.0

    # gradient of the objective against central differences
    cells = findall(vec(.!air_mask(start)))
    Gt = Sensitivity3D(start, d, PRISM; cells)
    o = Inv3DOptions()
    R = GravityGeophysics._inv3d_regularizer(start, cells, o)
    p = GravityGeophysics.Inv3DProblem(Gt, d.error, d.value ./ d.error, cells, R, sparse(R' * R),
                                       zeros(length(cells)), ones(length(cells)), o)
    m0 = 100 .* rand(MersenneTwister(3), length(cells))
    s = GravityGeophysics._inv3d_evaluate(p, m0, 3.0)
    g = GravityGeophysics.inv3d_gradient(p, s)
    objective(m) = GravityGeophysics._inv3d_evaluate(p, m, 3.0).objective
    for j in (1, 77, 300)
        e = zeros(length(cells)); e[j] = 1e-2
        @test g[j] ≈ (objective(m0 .+ e) - objective(m0 .- e)) / 2e-2 rtol = 1e-2
    end
end

@testset "3D INR inversion" begin
    truth, d = block_case()
    start = model_like(truth, zeros(size(truth.A)))

    config = INR3DConfig(width = 64, layers = 3, frequencies = 2, learning_rate = 1e-2, seed = 7)
    opt = Inv3DOptions(max_iter = 1500, target_rms = 1.0, bounds = (-500.0, 500.0), verbose = false)
    r = INR3D(start, d, PRISM; config, options = opt)
    @test r isa INR3DResult
    @test r.rms < 0.2 * r.history[1].rms
    @test all(-500 .<= r.density .<= 500)
    @test chi2_rms3d(d, data_like(d, r.predicted)).rms ≈ r.rms rtol = 1e-6
    @test r.beta == 0.0 && all(==(1.0), r.weights)

    # the returned model is the best-misfit network of the run
    @test r.rms ≈ minimum(h.rms for h in r.history) rtol = 1e-9
    @test r.reason != :target_rms || r.rms <= 1.0

    # the body is in the right place laterally; its depth is not asserted (no depth weighting)
    peak = Tuple(CartesianIndices(r.density)[argmax(r.density)])
    @test 3 <= peak[1] <= 7 && 4 <= peak[2] <= 8

    # a run stopped by the epoch budget reports the model it returns, not one Adam step behind it
    rs = INR3D(start, d, PRISM; config, options = Inv3DOptions(max_iter = 5, target_rms = 1.0,
                                                               bounds = (-500.0, 500.0), verbose = false))
    @test rs.reason == :max_iter && !rs.converged
    @test rs.history[end].iteration == 5
    @test rs.rms ≈ minimum(h.rms for h in rs.history) rtol = 1e-9

    # a diverging run stops and still returns the best network it saw
    rd = INR3D(start, d, PRISM;
               config = INR3DConfig(width = 64, layers = 3, learning_rate = 50.0, seed = 7),
               options = Inv3DOptions(max_iter = 200, target_rms = 1.0, bounds = (-500.0, 500.0),
                                      verbose = false))
    @test rd.reason == :diverged && !rd.converged
    @test rd.rms ≈ minimum(h.rms for h in rd.history) rtol = 1e-9
    @test rd.rms < rd.history[end].rms
end

# The rate schedule and the coarse-to-fine release of the Fourier bands, and switching them off
@testset "3D INR training schedule" begin
    alg = INR3DConfig(frequencies = 3, schedule = :cosine, rate_floor = 0.1, anneal = 0.5)
    rate(e) = GravityGeophysics._inr_rate(alg, 1.0, e, 1000)

    @test rate(20) == 1.0                                  # the peak lands at the end of the 2 % ramp
    @test rate(10) ≈ 0.5 && rate(1) ≈ 0.05                 # linear ramp before it
    @test rate(1000) ≈ 0.1 rtol = 1e-12                    # the floor at the last epoch
    @test rate(510) ≈ 0.55 atol = 0.01                     # half way down in the middle
    @test issorted([rate(e) for e in 20:1000]; rev = true) # and monotone in between
    @test all(e -> GravityGeophysics._inr_rate(INR3DConfig(schedule = :none), 1.0, e, 1000) == 1.0, 1:1000)

    # band k of the features is released over its own slice of the first `anneal` of the budget
    weights(e) = GravityGeophysics._inr_band_weights(alg, 3 + 6 * 3, e, 1000)
    @test weights(1)[1:3] == ones(3)                       # the raw coordinates always count
    @test weights(1)[4:9] ≈ fill(0.006, 6) atol = 1e-3     # the lowest band is just opening
    @test all(iszero, weights(1)[10:21])                   # the higher ones are shut
    @test weights(200)[4:9] == ones(6) && all(iszero, weights(200)[16:21])
    @test weights(500) == ones(21) && weights(1000) == ones(21)
    @test GravityGeophysics._inr_band_weights(INR3DConfig(frequencies = 3, anneal = 0.0), 21, 1, 1000) == ones(21)

    truth, d = block_case()
    start = model_like(truth, zeros(size(truth.A)))
    opt = Inv3DOptions(max_iter = 600, target_rms = 1.0, bounds = (-500.0, 500.0), verbose = false)
    base = (width = 64, layers = 3, frequencies = 4, learning_rate = 1e-2, seed = 7)
    r = INR3D(start, d, PRISM; config = INR3DConfig(; base...), options = opt)
    @test r.converged && r.rms <= 1.0
    peak = Tuple(CartesianIndices(r.density)[argmax(r.density)])
    @test 3 <= peak[1] <= 7 && 4 <= peak[2] <= 8

    # the switches are honoured: a constant rate with every band present is a different run
    rc = INR3D(start, d, PRISM; config = INR3DConfig(; base..., schedule = :none, anneal = 0.0), options = opt)
    @test rc.density != r.density

    @test_throws ArgumentError INR3D(start, d, PRISM; config = INR3DConfig(schedule = :linear), options = opt)
    @test_throws ArgumentError INR3D(start, d, PRISM; config = INR3DConfig(rate_floor = 0.0), options = opt)
    @test_throws ArgumentError INR3D(start, d, PRISM; config = INR3DConfig(anneal = 1.5), options = opt)
end

# Two blocks of opposite sign under noisy data, the one case whose target rms is the noise level
@testset "3D INR inversion -- twin bodies" begin
    truth, d = twin_case()
    start = model_like(truth, zeros(size(truth.A)))
    bounds = (-400.0, 400.0)

    config = INR3DConfig(width = 64, layers = 3, frequencies = 2, learning_rate = 1e-2, seed = 7)
    r = INR3D(start, d, PRISM; config,
              options = Inv3DOptions(; max_iter = 2000, target_rms = 1.0, bounds, verbose = false))

    # the network fits noisy data down to the noise level and stops there (discrepancy principle)
    @test r.converged && r.reason == :target_rms && r.rms <= 1.0
    @test r.history[1].rms > 10                                   # and had a long way to come
    @test length(r.history) - 1 < 2000                            # stopped on the misfit, not the budget
    @test chi2_rms3d(d, data_like(d, r.predicted)).rms ≈ r.rms rtol = 1e-6
    @test all(bounds[1] .<= r.density .<= bounds[2])

    # both bodies come back, with the right sign in the right corner and nothing much in the other two
    @test all(>(0), r.density[3:5, 3:5, 2:4])
    @test all(<(0), r.density[8:10, 8:10, 2:4])
    column = sum(r.density; dims = 3)[:, :, 1]                    # depth-integrated contrast, kg/m³·cell
    @test mean(column[3:5, 3:5]) > 200 && mean(column[8:10, 8:10]) < -200
    @test maximum(abs, column[3:5, 8:10]) < 100 && maximum(abs, column[8:10, 3:5]) < 100
    @test Tuple(CartesianIndices(r.density)[argmax(r.density)])[1:2] ⊆ 2:6
    @test Tuple(CartesianIndices(r.density)[argmin(r.density)])[1:2] ⊆ 7:11

    # as good a model as the Gauss-Newton run on the same data, which pays for depth weighting to get there
    g = Invert3D(start, d, PRISM;
                 options = Inv3DOptions(; max_iter = 30, target_rms = 1.0, bounds, weighting = :depth,
                                        verbose = false))
    @test g.reason == :target_rms
    @test cor(vec(r.density), vec(truth.A)) > 0.8 * cor(vec(g.density), vec(truth.A))

    # one network, one model: the same seed gives the same run back
    r2 = INR3D(start, d, PRISM; config,
               options = Inv3DOptions(; max_iter = 2000, target_rms = 1.0, bounds, verbose = false))
    @test r2.rms == r.rms && r2.density == r.density
end
