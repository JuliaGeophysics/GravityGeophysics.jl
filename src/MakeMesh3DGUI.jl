# Mesh designer window for MakeMesh3D
# Author: @pankajkmishra
# Map view of the mesh (cell lines, core box, stations) with the design parameters in text boxes, the
# checks under it, and buttons to update, zoom, save the start model and mask, and export the figure

#----- rectangle outline for the core box -----
corner_points(yl, yh, xl, xh) = [Point2f(yl, xl), Point2f(yh, xl), Point2f(yh, xh),
                                 Point2f(yl, xh), Point2f(yl, xl)]

const PARAM_HELP = Dict(
    "Core cell N–S (m)" =>
        "Core cell width in the N–S (X) direction, in metres. Use about half the " *
        "station spacing so no cell holds more than one station and structure between " *
        "stations is resolved. Smaller = finer resolution but more cells and slower.",
    "Core cell E–W (m)" =>
        "Core cell width in the E–W (Y) direction, in metres. Use about half the " *
        "station spacing (for line surveys, of the spacing along the lines).",
    "Core cells N–S" =>
        "Number of core cells in the N–S (X) direction. Independent of 'Core cell N–S (m)': " *
        "the core spans count × cell width, centred on the stations. An odd count keeps a " *
        "regular grid of stations on cell centres.",
    "Core cells E–W" =>
        "Number of core cells in the E–W (Y) direction. Independent of 'Core cell E–W (m)': " *
        "the core spans count × cell width, centred on the stations.",
    "Δρ background (kg/m³)" =>
        "Density contrast that fills the starting model, usually 0 (an anomaly inversion).",
    "Pad factor" =>
        "Geometric growth ratio of the padding cells (~1.3–1.5). The padding carries the " *
        "field of structure beyond the survey edge; potential fields decay as distance⁻², " *
        "so a padding as wide as the core is enough.",
    "Pad cells N–S" =>
        "Number of padding cells added on each side in N–S (X). The padding should be " *
        "at least as wide as the core.",
    "Pad cells E–W" =>
        "Number of padding cells added on each side in E–W (Y). The padding should be " *
        "at least as wide as the core.",
    "Layer (m)" =>
        "Thickness of the regular layers from the highest ground of the core down to the " *
        "depth of investigation; about half the station spacing (or the flight height).",
    "Vertical factor" =>
        "Geometric growth of the layers below the depth of investigation (~1.3–1.5), until " *
        "they are as thick as the regular block.",
    "Depth × core width" =>
        "Depth of investigation below the lowest ground of the core, as a multiple of the core " *
        "width. Gravity data resolve about half the survey width (~0.5).",
)

function _make_mesh3D_gui(;
    data_file, out_model, out_mask, topo_file, surface, lakes, uniform_layers,
    n_pad, pad_factor, vertical_factor, depth_mult,
    air_layers, air_factor, air_first_div,
    colormap, density_range, site_color, site_size_full, site_size_core,
    grid_color, grid_linewidth, fig_size, interactive,
    d, sx, sy, spacing, ρ_bg0, ρ_ref, dx_core0, dy_core0, nx_core0, ny_core0, z_first0, topo_ctx)

    use_topo = !isempty(strip(topo_file))
    pad0 = design_mesh(sx, sy; dx_core = dx_core0, dy_core = dy_core0, nx_core = nx_core0, ny_core = ny_core0,
                       nx_pad = n_pad, ny_pad = n_pad, pad_factor = pad_factor, z_first = z_first0,
                       z_factor = vertical_factor, depth_mult = depth_mult, uniform_layers = uniform_layers)

#----- draw the map view of the mesh -----
function draw2d!(ax, m)
    heatmap!(ax, m.y_edges_km, m.x_edges_km, fill(m.ρ_bg, m.ny, m.nx);
             colormap = colormap, colorrange = density_range)
    vlines!(ax, m.y_edges_km; color = grid_color, linewidth = grid_linewidth)
    hlines!(ax, m.x_edges_km; color = grid_color, linewidth = grid_linewidth)
    lines!(ax, corner_points(m.core_y0_km, m.core_y1_km, m.core_x0_km, m.core_x1_km);
           color = :black, linewidth = 2.5)
    scatter!(ax, sy ./ 1000, sx ./ 1000; color = site_color, marker = :circle,
             markersize = site_size)
    return ax
end

fig = Figure(size = fig_size, figure_padding = (40, 14, 14, 14))

status = Observable("Outputs → $(basename(out_model)), $(basename(out_mask))")
help_obs = Observable("Click  ?  next to a parameter for an explanation.")
sug_obs = Observable("")
sug_color = Observable(RGBf(0.15, 0.5, 0.25))
show_full = Observable(true)
site_size = Observable(site_size_full)

Colorbar(fig[1, 3]; colormap = colormap, colorrange = density_range,
    label = "Δρ (kg/m³)", width = 16)

left = GridLayout(fig[1, 1]; valign = :top)

infotext = join([
    @sprintf("datafile  >  %s", basename(data_file)),
    @sprintf("stations  >  %d", length(sx)),
    @sprintf("data      >  %s (%s)", _COMPONENT_NAME[d.component], units(d.component)),
    @sprintf("spacing   >  %.0f m", spacing.median),
    @sprintf("topo      >  %s", use_topo ? basename(topo_file) : surface !== nothing ? "surface function" : "off"),
    @sprintf("mask(out) >  %s", basename(out_mask)),
], "\n")
Label(left[1, 1:3], infotext; font = "Consolas", fontsize = 13, halign = :left,
    justification = :left, color = :gray25, tellwidth = false)

next_row = Ref(1)

function add_param!(name, default; isint = false)
    r = (next_row[] += 1)
    Label(left[r, 1], name; halign = :right, fontsize = 13)
    s = isint ? string(Int(round(default))) : string(default)
    tb = Textbox(left[r, 2]; stored_string = s, validator = Float64, width = 84)
    ib = Button(left[r, 3]; label = "?", fontsize = 13, width = 26)
    on(ib.clicks) do _
        help_obs[] = get(PARAM_HELP, name, "")
    end
    return tb
end

tb_dx   = add_param!("Core cell N–S (m)", dx_core0)
tb_dy   = add_param!("Core cell E–W (m)", dy_core0)
tb_nxc  = add_param!("Core cells N–S", nx_core0; isint = true)
tb_nyc  = add_param!("Core cells E–W", ny_core0; isint = true)
tb_kap  = add_param!("Δρ background (kg/m³)", ρ_bg0)
tb_npx  = add_param!("Pad cells N–S", pad0.nx_pad; isint = true)
tb_npy  = add_param!("Pad cells E–W", pad0.ny_pad; isint = true)
tb_pf   = add_param!("Pad factor", pad_factor)
tb_zf   = add_param!("Layer (m)", z_first0)
tb_zfac = add_param!("Vertical factor", vertical_factor)
tb_dm   = add_param!("Depth × core width", depth_mult)

updatebtn = Button(left[(next_row[] += 1), 1:3]; label = "Update mesh", fontsize = 14)

viewgrid = GridLayout(left[(next_row[] += 1), 1:3])
corefullbtn = Button(viewgrid[1, 1]; label = "Show core", fontsize = 14, tellwidth = false)
resetbtn    = Button(viewgrid[1, 2]; label = "Reset zoom", fontsize = 14, tellwidth = false)
colsize!(viewgrid, 1, Relative(0.5)); colsize!(viewgrid, 2, Relative(0.5))
colgap!(viewgrid, 8)

actiongrid = GridLayout(left[(next_row[] += 1), 1:3])
savebtn   = Button(actiongrid[1, 1]; label = "Save model + mask", fontsize = 14, tellwidth = false)
exportbtn = Button(actiongrid[1, 2]; label = "Export figure", fontsize = 14, tellwidth = false)
colsize!(actiongrid, 1, Relative(0.5)); colsize!(actiongrid, 2, Relative(0.5))
colgap!(actiongrid, 8)

rowgap!(left, 7)
colgap!(left, 8)
colsize!(left, 2, Fixed(86))
colsize!(left, 3, Fixed(26))

Label(fig[2, 1:3], sug_obs; fontsize = 13, halign = :left, justification = :left,
    color = sug_color, word_wrap = true, tellwidth = false)
Label(fig[3, 1:3], help_obs; fontsize = 12.5, halign = :left, justification = :left,
    color = :gray30, word_wrap = true, tellwidth = false)
Label(fig[4, 1:3], status; fontsize = 11, halign = :left, color = :gray45,
    word_wrap = true, tellwidth = false)

colsize!(fig.layout, 1, Fixed(290))
colsize!(fig.layout, 3, Fixed(70))
rowsize!(fig.layout, 1, Relative(0.82))
for r in 2:4
    rowsize!(fig.layout, r, Auto())
end
rowgap!(fig.layout, 6)

tbget(tb, dflt) = begin
    s = tb.displayed_string[]
    (s === nothing || isempty(strip(s))) && return dflt
    v = tryparse(Float64, s)
    v === nothing ? dflt : v
end

function build_mesh()
    build_mesh_bundle(d, sx, sy;
        ρ_bg = tbget(tb_kap, ρ_bg0), ρ_ref = ρ_ref,
        dx_core = tbget(tb_dx, dx_core0),
        dy_core = tbget(tb_dy, dy_core0),
        nx_core = Int(round(tbget(tb_nxc, Float64(nx_core0)))),
        ny_core = Int(round(tbget(tb_nyc, Float64(ny_core0)))),
        nx_pad = Int(round(tbget(tb_npx, Float64(pad0.nx_pad)))),
        ny_pad = Int(round(tbget(tb_npy, Float64(pad0.ny_pad)))),
        pad_factor = tbget(tb_pf, pad_factor),
        z_first = tbget(tb_zf, z_first0),
        z_factor = tbget(tb_zfac, vertical_factor),
        depth_mult = tbget(tb_dm, depth_mult), uniform_layers = uniform_layers,
        topo_ctx = topo_ctx, topo_path = topo_file, surface = surface, lakes = lakes,
        n_air = air_layers, air_factor = air_factor, air_first_div = air_first_div)
end

mesh = Observable(build_mesh())
current_ax = Ref{Any}(nothing)

function apply_limits!(ax, m)
    if show_full[]
        xlims!(ax, m.full_y0_km, m.full_y1_km)
        ylims!(ax, m.full_x0_km, m.full_x1_km)
    else
        padx = 0.05 * (m.core_y1_km - m.core_y0_km)
        pady = 0.05 * (m.core_x1_km - m.core_x0_km)
        xlims!(ax, m.core_y0_km - padx, m.core_y1_km + padx)
        ylims!(ax, m.core_x0_km - pady, m.core_x1_km + pady)
    end
    return nothing
end

function refresh!(m)
    current_ax[] === nothing || delete!(current_ax[])
    ax = Axis(fig[1, 2]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right,
        xgridvisible = false, ygridvisible = false,
        xlabelfont = :bold, ylabelfont = :bold,
        xticklabelcolor = :gray55, yticklabelcolor = :gray55)
    draw2d!(ax, m)
    current_ax[] = ax
    apply_limits!(ax, m)
    s = suggestions(m)
    topo_info = m.topo_active ? @sprintf("topo %s · air %d", m.topo_name, m.nz_air) : "topo flat"
    site_info = m.site_ground.ok ?
        (m.site_surface.nadjusted == 0 ? "all stations above ground" : @sprintf("all stations above ground (%d lowered column(s))", m.site_surface.nadjusted)) :
        @sprintf("%d station(s) in the ground", m.site_ground.nbad)
    gridinfo = @sprintf("%d × %d × %d cells · core %d × %d (%.0f × %.0f m) · max stations/cell %d · %d regular layers · depth %.2f km (investigation %.2f km) · %s · %s",
        m.nx, m.ny, m.nz, m.nx_core, m.ny_core, m.dx_core, m.dy_core, m.maxper, m.nz_regular, m.depth_km, m.doi_km, topo_info, site_info)
    sug_obs[] = gridinfo * "\n" * s.text
    sug_color[] = s.ok ? RGBf(0.15, 0.5, 0.25) : RGBf(0.75, 0.2, 0.15)
    return nothing
end

on(mesh) do m
    refresh!(m)
end
refresh!(mesh[])

on(updatebtn.clicks) do _
    mesh[] = build_mesh()
end

for tb in (tb_dx, tb_dy, tb_nxc, tb_nyc, tb_kap, tb_npx, tb_npy, tb_pf, tb_zf, tb_zfac, tb_dm)
    on(tb.stored_string) do _
        mesh[] = build_mesh()
    end
end

on(corefullbtn.clicks) do _
    show_full[] = !show_full[]
    corefullbtn.label[] = show_full[] ? "Show core" : "Show full"
    site_size[] = show_full[] ? site_size_full : site_size_core
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(resetbtn.clicks) do _
    current_ax[] === nothing || apply_limits!(current_ax[], mesh[])
end

on(savebtn.clicks) do _
    m = mesh[]
    _save_mesh_outputs(m, out_model, out_mask)
    ground_status = m.site_ground.ok ?
        (m.site_surface.nadjusted == 0 ? "all stations above ground" : @sprintf("all stations above ground after %d lowered column(s)", m.site_surface.nadjusted)) :
        @sprintf("WARNING: %d station(s) in the ground; see guidance above", m.site_ground.nbad)
    status[] = @sprintf("Saved %d×%d×%d model → %s | mask → %s | %s", m.nx, m.ny, m.nz, out_model, out_mask, ground_status)
    @info status[]
end

on(exportbtn.clicks) do _
    m = mesh[]
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    efig = Figure(size = (1200, 1000))
    eax = Axis(efig[1, 1]; xlabel = "Y East (km)", ylabel = "X North (km)",
        aspect = DataAspect(), halign = :right, xgridvisible = false, ygridvisible = false)
    draw2d!(eax, m)
    apply_limits!(eax, m)
    Colorbar(efig[1, 2]; colormap = colormap, colorrange = density_range,
        label = "Δρ (kg/m³)", width = 16)
    fn = joinpath(dirname(out_model), "mesh_map_$ts.png")
    save(fn, efig; px_per_unit = 3)
    status[] = "Saved figure → $fn"
    @info status[]
end

    interactive && _open_window(fig)

    return fig, mesh
end
