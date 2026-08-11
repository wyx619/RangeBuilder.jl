struct PlotCommand
    kind::Symbol
    data::Any
    color::Any
    linewidth::Float64
    linestyle::Symbol
    text::Union{Nothing,Vector{String}}
end

struct PlotSpec
    commands::Vector{PlotCommand}
    xlim::NTuple{2,Float64}
    ylim::NTuple{2,Float64}
end

function _plot_mode(wlines)
    mode = wlines isa Symbol ? wlines : Symbol(wlines)
    mode in (:none, :both, :del, :vor) ||
        throw(ArgumentError("wlines must be :none, :both, :del, or :vor"))
    return mode
end

function _plot_values(value, n, default)
    values = value === nothing ? default : (value isa AbstractVector ? collect(value) : [value])
    isempty(values) && throw(ArgumentError("plot style values must not be empty"))
    return [values[mod1(i, length(values))] for i in 1:n]
end

function _plot_limits(x, xlim, ylim)
    xr = xlim === nothing ? (minimum(x[:, 1]), maximum(x[:, 1])) : (Float64(xlim[1]), Float64(xlim[2]))
    yr = ylim === nothing ? (minimum(x[:, 2]), maximum(x[:, 2])) : (Float64(ylim[1]), Float64(ylim[2]))
    return xr, yr
end

function _r_colors(n)
    palette = (:black, :red, :green, :blue, :cyan, :magenta)
    return [palette[mod1(i, length(palette))] for i in 1:n]
end

function _plot_labels(x, color)
    dx = 0.02 * (maximum(x[:, 1]) - minimum(x[:, 1]))
    dy = 0.02 * (maximum(x[:, 2]) - minimum(x[:, 2]))
    positions = copy(x)
    positions[:, 1] .+= dx
    positions[:, 2] .+= dy
    return PlotCommand(:labels, positions, color, 0.0, :solid, string.(axes(x, 1)))
end

function _delvor_commands(dv, mode, points_color, del_color, vor_color, label_color,
                          linewidth, wpoints, number)
    commands = PlotCommand[]
    x = dv.x
    wpoints && push!(commands, PlotCommand(:points, x, points_color, linewidth, :solid, nothing))
    number && push!(commands, _plot_labels(x, label_color))
    if mode in (:both, :del)
        for i in axes(dv.mesh, 1)
            push!(commands, PlotCommand(:segment, [dv.mesh[i, 3] dv.mesh[i, 4];
                                                    dv.mesh[i, 5] dv.mesh[i, 6]],
                                        del_color, linewidth, :solid, nothing))
        end
    end
    if mode in (:both, :vor)
        for i in axes(dv.mesh, 1)
            linestyle = (dv.mesh[i, 11] == 1 || dv.mesh[i, 12] == 1) ? :dash : :solid
            push!(commands, PlotCommand(:segment, [dv.mesh[i, 7] dv.mesh[i, 8];
                                                    dv.mesh[i, 9] dv.mesh[i, 10]],
                                        vor_color, linewidth, linestyle, nothing))
        end
    end
    return commands
end

function _ashape_plotdata(sh::AShape; wlines=:none, wpoints=true, number=false,
                          col=nothing, lwd=nothing, xlim=nothing, ylim=nothing)
    mode = _plot_mode(wlines)
    colors = _plot_values(col, 5, _r_colors(5))
    widths = Float64.(_plot_values(lwd, 2, (1.0, 2.0)))
    commands = PlotCommand[]
    if mode == :none
        x = sh.x
        wpoints && push!(commands, PlotCommand(:points, x, colors[2], widths[1], :solid, nothing))
        number && push!(commands, _plot_labels(x, colors[5]))
    else
        append!(commands, _delvor_commands(sh.delvor, mode, colors[2], colors[3], colors[4],
                                           colors[5], widths[1], wpoints, number))
    end
    for i in axes(sh.edges, 1)
        push!(commands, PlotCommand(:segment, [sh.edges[i, 3] sh.edges[i, 4];
                                                sh.edges[i, 5] sh.edges[i, 6]],
                                    colors[1], widths[2], :solid, nothing))
    end
    xr, yr = _plot_limits(sh.x, xlim, ylim)
    return PlotSpec(commands, xr, yr)
end

function _ahull_plotdata(h::AHull; do_shape=false, wlines=:none, wpoints=true,
                         number=false, col=nothing, lwd=nothing, xlim=nothing, ylim=nothing)
    mode = _plot_mode(wlines)
    colors = _plot_values(col, 6, _r_colors(6))
    widths = Float64.(_plot_values(lwd, 3, (1.0, 1.0, 2.0)))
    if do_shape
        spec = _ashape_plotdata(h.ashape; wlines=mode, wpoints=wpoints, number=number,
                                col=colors[2:6], lwd=widths[1:2], xlim=xlim, ylim=ylim)
        commands = copy(spec.commands)
    elseif mode == :none
        commands = PlotCommand[]
        x = h.ashape.x
        wpoints && push!(commands, PlotCommand(:points, x, colors[3], widths[1], :solid, nothing))
        number && push!(commands, _plot_labels(x, colors[6]))
    else
        commands = _delvor_commands(h.ashape.delvor, mode, colors[3], colors[4], colors[5],
                                    colors[6], widths[1], wpoints, number)
    end
    arc_rows = findall(>(0), h.arcs[:, 3])
    for i in arc_rows
        push!(commands, PlotCommand(:path,
                                    arc(h.arcs[i, 1:2], h.arcs[i, 3], h.arcs[i, 4:5], h.arcs[i, 6]),
                                    colors[1], widths[3], :solid, nothing))
    end
    point_rows = findall(==(0), h.arcs[:, 3])
    if !isempty(point_rows)
        data = h.arcs[point_rows, 1:2]
        push!(commands, PlotCommand(:points, data, colors[1], widths[3], :solid, nothing))
    end
    xr, yr = _plot_limits(h.ashape.x, xlim, ylim)
    return PlotSpec(commands, xr, yr)
end

"""Build backend-independent drawing commands with R `plot.*` semantics."""
plotdata(dv::DelVor; wlines=:both, wpoints=true, number=false, col=nothing,
         lwd=1.0, xlim=nothing, ylim=nothing) = begin
    mode = _plot_mode(wlines)
    colors = _plot_values(col, 4, _r_colors(4))
    width = Float64(_plot_values(lwd, 1, (1.0,))[1])
    PlotSpec(_delvor_commands(dv, mode, colors[1], colors[2], colors[3], colors[4],
                              width, wpoints, number), _plot_limits(dv.x, xlim, ylim)...)
end

plotdata(sh::AShape; kwargs...) = _ashape_plotdata(sh; kwargs...)
plotdata(h::AHull; kwargs...) = _ahull_plotdata(h; kwargs...)
