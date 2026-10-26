# Sensitivity matrix
# Author: @pankajkmishra
# Gᵀ against the forward solver, the active-cell subset, the memory-mapped variant and the file workflow of
# ForwardSolve3D

using Test
using LinearAlgebra

@testset "sensitivity" begin
    mktempdir() do dir
        m = small_mesh(n = 6, nz = 4, h = 40.0)
        m.A[3:4, 3:4, 2:3] .= 150.0
        d = small_survey(n = 5, h = 35.0, height = 15.0)

        Gt = Sensitivity3D(m, d)
        @test size(Gt) == (length(m.A), ndata(d)) && eltype(Gt) == Float32
        @test Float64.(Gt' * Float32.(vec(m.A))) ≈ forward3d(m, d) rtol = 1e-4

        # a subset of cells gives the field of those cells alone
        cells = findall(vec(m.A) .!= 0)
        Gc = Sensitivity3D(m, d; cells)
        @test size(Gc) == (length(cells), ndata(d))
        @test Gc ≈ Gt[cells, :]

        # the sensitivity is the derivative of the response: doubling a cell doubles its contribution
        m2 = model_like(m, copy(m.A))
        m2.A[cells[1]] += 100.0
        @test forward3d(m2, d) ≈ forward3d(m, d) .+ 100 .* Float64.(Gt[cells[1], :]) rtol = 1e-5

        # with air cells the response is that of the earth cells alone
        topo = model_like(m, copy(m.A))
        topo.A[:, :, 1] .= NaN
        Gair = Sensitivity3D(topo, d; cells = findall(vec(.!air_mask(topo))))
        @test Float64.(Gair' * Float32.(filter(!isnan, vec(topo.A)))) ≈ forward3d(topo, d) rtol = 1e-4

        # saved and reloaded as jld2, the matrix and its cell indices come back unchanged
        gpath = SaveSensitivity3D(joinpath(dir, "sensitivity.jld2"), Gt)
        back = LoadSensitivity3D(gpath)
        @test back.sensitivity ≈ Gt && back.cells == collect(1:length(m.A))
        @test_throws ErrorException LoadSensitivity3D(joinpath(dir, "missing.jld2"))

        # ForwardSolve3D writes the predicted data of a model file onto the survey of a data file
        mpath = WriteModel3D(joinpath(dir, "model.den"), m)
        dpath = WriteData3D(joinpath(dir, "data.xyz"), d)
        fpath = WriteFwdCtrl3D(joinpath(dir, "FwdCtrl"), FwdCtrl3D())
        pred = ReadData3D(ForwardSolve3D(mpath, dpath, fpath))
        @test pred.value ≈ forward3d(m, d) rtol = 1e-4
        @test pred.x == d.x && pred.error == d.error
    end
end
