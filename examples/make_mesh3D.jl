using GravityGeophysics

data_file = length(ARGS) >= 1 ? ARGS[1] : ""
out_model = length(ARGS) >= 2 ? ARGS[2] : (isempty(data_file) ? "" : joinpath(dirname(abspath(data_file)), "mesh_start_model.den"))
topo_file = length(ARGS) >= 3 ? ARGS[3] : ""
out_mask  = length(ARGS) >= 4 ? ARGS[4] : (isempty(out_model) ? "" : joinpath(dirname(out_model), "mesh_start_model.mask"))
ρ_bg      = length(ARGS) >= 5 ? something(tryparse(Float64, ARGS[5]), 0.0) : 0.0
mode      = length(ARGS) >= 6 ? Symbol(lowercase(ARGS[6])) : Symbol(lowercase(get(ENV, "MAGGEO_MESH_MODE", "gui")))

topo_crs        = "EPSG:4326"
cell_width_frac = 0.5
n_pad           = nothing
pad_factor      = 1.4
first_layer_div = 2.0
vertical_factor = 1.4
depth_mult      = 0.5
air_layers      = 0
air_factor      = 1.35
air_first_div   = 2.0

colormap             = COLORMAP
density_range = (-300.0, 300.0)
site_color           = :black
site_size_full       = 6
site_size_core       = 5
grid_color           = (:grey, 0.7)
grid_linewidth       = 1.0
fig_size             = (1500, 1000)

MakeMesh3D(data_file;
    out_model            = out_model,
    out_mask             = out_mask,
    topo_file            = topo_file,
    topo_crs             = topo_crs,
    ρ_bg                 = ρ_bg,
    mode                 = mode,
    cell_width_frac      = cell_width_frac,
    n_pad                = n_pad,
    pad_factor           = pad_factor,
    first_layer_div      = first_layer_div,
    vertical_factor      = vertical_factor,
    depth_mult           = depth_mult,
    air_layers           = air_layers,
    air_factor           = air_factor,
    air_first_div        = air_first_div,
    colormap             = colormap,
    density_range = density_range,
    site_color           = site_color,
    site_size_full       = site_size_full,
    site_size_core       = site_size_core,
    grid_color           = grid_color,
    grid_linewidth       = grid_linewidth,
    fig_size             = fig_size,
)
