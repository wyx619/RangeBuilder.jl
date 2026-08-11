module RangeBuilderMakieExt

using RangeBuilder
import Makie

function _draw!(ax, spec::RangeBuilder.PlotSpec)
    for command in spec.commands
        if command.kind == :points
            data = command.data
            Makie.scatter!(ax, data[:, 1], data[:, 2], color=command.color, markersize=8)
        elseif command.kind == :segment || command.kind == :path
            data = command.data
            Makie.lines!(ax, data[:, 1], data[:, 2], color=command.color,
                         linewidth=command.linewidth, linestyle=command.linestyle)
        elseif command.kind == :labels
            data = command.data
            Makie.text!(ax, data[:, 1], data[:, 2], text=command.text,
                        color=command.color, align=(:left, :bottom))
        end
    end
    Makie.xlims!(ax, spec.xlim...)
    Makie.ylims!(ax, spec.ylim...)
    return ax
end

function Makie.plot!(ax::Makie.Axis, dv::RangeBuilder.DelVor; kwargs...)
    _draw!(ax, RangeBuilder.plotdata(dv; kwargs...))
end

function Makie.plot!(ax::Makie.Axis, sh::RangeBuilder.AShape; kwargs...)
    _draw!(ax, RangeBuilder.plotdata(sh; kwargs...))
end

function Makie.plot!(ax::Makie.Axis, h::RangeBuilder.AHull; kwargs...)
    _draw!(ax, RangeBuilder.plotdata(h; kwargs...))
end

function Makie.plot(dv::RangeBuilder.DelVor; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
    Makie.plot!(ax, dv; kwargs...)
    return fig
end

function Makie.plot(sh::RangeBuilder.AShape; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
    Makie.plot!(ax, sh; kwargs...)
    return fig
end

function Makie.plot(h::RangeBuilder.AHull; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
    Makie.plot!(ax, h; kwargs...)
    return fig
end

end
