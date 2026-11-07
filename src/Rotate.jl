# Data rotation
# Author: @pankajkmishra
# Turn the data frame by θ degrees clockwise, for a rotated mesh or a strike direction, and keep the
# history in the file. The mesh is never rotated (the prism solution needs axis-aligned cells).
# Only station positions and the frame azimuth change: gz does not depend on the frame.

const ROTATION_KINDS = (:mesh, :strike)

"""
    rotate_data(d, θ; kind=:mesh) -> Data

Turn the frame of `d` by θ degrees clockwise: x' = x cos θ + y sin θ, y' = −x sin θ + y cos θ.
`kind` (`:mesh` or `:strike`) is kept in the rotation history. Rotate by −θ to undo.
"""
function rotate_data(d::Data, θ::Real; kind::Symbol = :mesh)
    kind in ROTATION_KINDS || throw(ArgumentError("rotation kind must be one of $ROTATION_KINDS"))
    c, s = cosd(θ), sind(θ)
    out = deepcopy(d)
    out.x, out.y = c .* d.x .+ s .* d.y, -s .* d.x .+ c .* d.y
    out.azimuth = rem(d.azimuth + θ, 360, RoundNearest)
    out.rotations = vcat(d.rotations, (kind, Float64(θ)))
    out
end
