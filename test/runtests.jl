# Package tests
# Author: @pankajkmishra
# Runs every test file: model, mask and data files, forward modelling and its accuracy, the deterministic and
# inversions, viewers, mesh design, editors and mesh-to-mesh projection
# Pkg.test(test_args = ["accuracy"]) runs the accuracy checks only

using Test
using GravityGeophysics
using Printf
using LinearAlgebra
using Random
using Statistics
using CairoMakie: set_close_to!

const DATA = joinpath(@__DIR__, "data")
include(joinpath(@__DIR__, "Helpers3D.jl"))

@testset "GravityGeophysics.jl Tests" begin
    if "accuracy" in ARGS
        include(joinpath(@__DIR__, "TestAccuracy3D.jl"))
    else
        #---------- files ----------
        include(joinpath(@__DIR__, "TestIO3D.jl"))
        include(joinpath(@__DIR__, "TestMask3D.jl"))
        include(joinpath(@__DIR__, "TestData3D.jl"))

        #---------- forward ----------
        include(joinpath(@__DIR__, "TestFwd3D.jl"))
        include(joinpath(@__DIR__, "TestAccuracy3D.jl"))

        #---------- inversion ----------
        include(joinpath(@__DIR__, "TestInv3D.jl"))

        #---------- viewers, mesh design and model tools ----------
        include(joinpath(@__DIR__, "TestPlot3D.jl"))
        include(joinpath(@__DIR__, "TestMesh3D.jl"))
        include(joinpath(@__DIR__, "TestEdit3D.jl"))
        include(joinpath(@__DIR__, "TestMeshToMesh3D.jl"))
    end
end
