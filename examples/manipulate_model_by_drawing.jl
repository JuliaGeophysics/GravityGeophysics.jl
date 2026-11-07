using GravityGeophysics

model_file = length(ARGS) >= 1 ? ARGS[1] : ""

replacement_density = 300.0
layers_above               = 2
layers_below               = 2
transition_layers          = 1
apply_to_all_depths        = false
depth_range                = (0.0, 50000.0)

log10_scale    = false
colormap       = COLORMAP
with_padding   = true
max_depth      = nothing
value_range    = (0.0, 300.0)
show_grid      = true
grid_color     = :black
grid_linewidth = 0.5
grid_alpha     = 0.3
pad_tol        = 0.5

viewer_figsize = (1200, 950)

EditModelByDrawing(model_file;
    replacement_density = replacement_density,
    layers_above               = layers_above,
    layers_below               = layers_below,
    transition_layers          = transition_layers,
    apply_to_all_depths        = apply_to_all_depths,
    depth_range                = depth_range,
    log10_scale                = log10_scale,
    colormap                   = colormap,
    with_padding               = with_padding,
    max_depth                  = max_depth,
    value_range                = value_range,
    show_grid                  = show_grid,
    grid_color                 = grid_color,
    grid_linewidth             = grid_linewidth,
    grid_alpha                 = grid_alpha,
    pad_tol                    = pad_tol,
    viewer_figsize             = viewer_figsize,
)
