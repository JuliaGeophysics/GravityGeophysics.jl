# Author: @pankajkmishra

using Test
using GravityGeophysics

include("Helpers3D.jl")

@testset "GravityGeophysics.jl" begin
    include("TestIO3D.jl")
    include("TestControl3D.jl")
    include("TestForward3D.jl")
    include("TestSensitivity3D.jl")
    include("TestMesh3D.jl")
    include("TestInversion3D.jl")
    include("TestWorkflows3D.jl")
end
