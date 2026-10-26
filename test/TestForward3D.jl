# 3D gravity forward solvers
# Author: @pankajkmishra
# The prism sum against the analytic field of a uniform sphere (a sphere attracts like a point mass outside
# itself), linearity in density, the finite-difference Poisson solver against the prism sum, and the receiver
# check that keeps observation points out of earth cells

using Test
using LinearAlgebra
using Statistics

@testset "3D forward solvers" begin
    # a sphere discretised on a fine mesh, its field read 100 m above the centre and off to the sides
    h, n, nz = 20.0, 24, 12
    m = GravModel3D(fill(h, n), fill(h, n), fill(h, nz), zeros(n, n, nz), [-n * h / 2, -n * h / 2, 0.0])
    centre, radius, drho = (0.0, 0.0, 120.0), 80.0, 400.0
    for k in eachindex(m.cz), j in eachindex(m.cy), i in eachindex(m.cx)
        (m.cx[i] - centre[1])^2 + (m.cy[j] - centre[2])^2 + (m.cz[k] - centre[3])^2 <= radius^2 &&
            (m.A[i, j, k] = drho)
    end
    d = small_survey(n = 7, h = 60.0, height = 40.0)
    exact = sphere_anomaly(d, centre, radius, drho)
    prism = forward3d(m, d)
    # the staircase sphere holds about the right mass, so the fields agree to a few per cent
    @test maximum(abs, prism .- exact) / maximum(abs, exact) < 0.05
    @test cor(prism, exact) > 0.999

    # the solver is linear in the density contrast and zero for an empty model
    @test forward3d(model_like(m, 2 .* m.A), d) ≈ 2 .* prism
    @test all(==(0), forward3d(model_like(m, zeros(size(m.A))), d))

    # a sphere above the sensors pulls the opposite way (gz is positive down)
    below = small_survey(n = 7, h = 60.0, height = 40.0)
    below.elevation .= -400.0
    @test all(<(0), forward3d(m, below))

    # finite differences agree with the prism sum on a uniform mesh, with the sensors in the air above it
    relrms(a, b) = sqrt(mean(abs2, a .- b)) / sqrt(mean(abs2, b))
    fd = forward3d_fd(m, d, FDCtrl3D(tol = 1e-11, maxiter = 20_000))
    @test all(isfinite, fd)
    @test relrms(fd, prism) < 0.05
    @test cor(fd, prism) > 0.99

    # air cells under topography carry no mass and the sensors above them stay inside the padded domain
    tm = model_like(m, copy(m.A))
    tm.A[:, :, 1] .= NaN
    tm.A[1:(n ÷ 2), :, 2] .= NaN
    @test relrms(forward3d_fd(tm, d), forward3d(tm, d)) < 0.05

    # a stretched, laterally padded mesh is accepted as it is, the case MakeMesh3D writes
    w = vcat(reverse([h * 1.5^i for i in 1:4]), fill(h, n), [h * 1.5^i for i in 1:4])
    v = vcat(fill(h, nz), [h * 1.3^i for i in 1:4])
    sm = GravModel3D(w, w, v, zeros(length(w), length(w), length(v)),
                     [-sum(w) / 2, -sum(w) / 2, 0.0])
    for k in eachindex(sm.cz), j in eachindex(sm.cy), i in eachindex(sm.cx)
        (sm.cx[i] - centre[1])^2 + (sm.cy[j] - centre[2])^2 + (sm.cz[k] - centre[3])^2 <= radius^2 &&
            (sm.A[i, j, k] = drho)
    end
    @test relrms(forward3d_fd(sm, d), forward3d(sm, d)) < 0.10

    # a sensor inside an earth cell is refused, one in an air cell is fine
    inside = small_survey(n = 2, h = 50.0)
    inside.elevation .= -100.0
    @test_throws ErrorException forward3d(m, inside)
    topo = model_like(m, copy(m.A))
    topo.A[:, :, 1] .= NaN
    air_level = small_survey(n = 2, h = 50.0)
    air_level.elevation .= -10.0
    @test check_receivers(topo, air_level) === nothing
end
