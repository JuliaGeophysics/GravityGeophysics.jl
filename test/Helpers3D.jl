# Shared 3D test fixtures
# Author: @pankajkmishra
# A small survey and mesh used by several test files
# the sensor height is an elevation above the mesh top (depth zero), not a negative z

# n × n grid of sensors, spacing h, elevation above the top of the mesh
function small_survey(; n = 9, h = 50.0, height = 20.0)
    c = ((1:n) .- (n + 1) / 2) .* h
    GravData3D(x = repeat(c, inner = n), y = repeat(c, outer = n), elevation = fill(height, n * n),
               value = zeros(n * n), error = ones(n * n))
end

# uniform cubes of size h, top at elevation 0, centred on the origin
small_mesh(; n = 12, nz = 6, h = 50.0, rho = 0.0) =
    GravModel3D(fill(h, n), fill(h, n), fill(h, nz), fill(rho, n, n, nz), [-n * h / 2, -n * h / 2, 0.0])

# a buried block of `rho` in the middle of a `small_mesh`, and the data it produces
function block_case(; n = 10, nz = 6, h = 50.0, rho = 300.0, nd = 11, hd = 40.0, floor = 0.01, percent = 2.0)
    truth = small_mesh(; n, nz, h)
    truth.A[4:6, 5:7, 2:4] .= rho
    d = small_survey(n = nd, h = hd)
    d.value .= forward3d(truth, d)
    d.error .= floor .+ percent / 100 .* abs.(d.value)
    truth, d
end

# two blocks of opposite sign in a `small_mesh` and their noisy data, the noise written as the error
function twin_case(; n = 12, nz = 6, h = 50.0, rho = 300.0, nd = 13, hd = 50.0,
                     floor = 0.005, percent = 2.0, seed = 11)
    truth = small_mesh(; n, nz, h)
    truth.A[3:5, 3:5, 2:4] .= rho
    truth.A[n-4:n-2, n-4:n-2, 2:4] .= -rho
    d = small_survey(n = nd, h = hd)
    d.value .= forward3d(truth, d)
    truth, add_noise(d; floor, percent, seed)
end

const PRISM = FwdCtrl3D()
