using Test
using RangeBuilder
using CairoMakie

points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
delvor_object = delvor(points)
shape = ashape(delvor_object; alpha=0.7)
hull = ahull(delvor_object; alpha=0.7)

@testset "Makie extension" begin
    @test !isnothing(Base.get_extension(RangeBuilder, :RangeBuilderMakieExt))

    figure = CairoMakie.plot(hull; do_shape=true, wlines=:vor, number=true)
    @test figure isa CairoMakie.Figure

    output = joinpath(mktempdir(), "ahull.png")
    CairoMakie.save(output, figure)
    @test isfile(output)
    @test filesize(output) > 0

    figure = CairoMakie.Figure()
    axis = CairoMakie.Axis(figure[1, 1])
    @test CairoMakie.plot!(axis, delvor_object; wlines=:both) === axis
    @test CairoMakie.plot!(axis, shape; wlines=:none) === axis
end
