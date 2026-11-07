# Shared test helpers
# Author: @pankajkmishra
# Station data in frame coordinates and the report line of a measured error against its tolerance

const GG = GravityGeophysics

function report(name, err, tol)
    @printf("  %-58s %10.2e  (tol %.0e)\n", name, err, tol)
    @test err <= tol
end

# gz stations in frame coordinates (x north, y east, z down)
station_data(x, y, z) =
    Data(collect(Float64, x), collect(Float64, y), collect(Float64, z), fill(NaN, length(x)), fill(NaN, length(x)),
         zeros(length(x)), ones(length(x)), string.(eachindex(x)), :gz, [NaN, NaN], 0.0, 0.0,
         Tuple{Symbol,Float64}[], "frame coordinates", 0.0, "test")
