# Optional differential verification for the first SHull porting stage.
# Run explicitly with: julia --project=test test/shull_prototype.jl
using Random
using RCall
using Test
using RangeBuilder

function r_shull_edges(points)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y
    RCall.reval("library(interp)")
    mesh = RCall.rcopy(RCall.reval("tri.mesh(x, y)\$arcs"))
    return Set((min(Int(row[1]), Int(row[2])), max(Int(row[1]), Int(row[2])))
               for row in eachrow(mesh))
end

function r_shull_hull(points)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y
    RCall.reval("library(interp)")
    return Set(Int.(RCall.rcopy(RCall.reval("tri.mesh(x, y)\$chull"))))
end

function r_shull_trlist(points)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y
    RCall.reval("library(interp)")
    return Int.(RCall.rcopy(RCall.reval("tri.mesh(x, y)\$trlist")))
end

function r_delvor_mesh(points)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y
    RCall.reval("library(alphahull)")
    return Float64.(RCall.rcopy(RCall.reval("alphahull::delvor(cbind(x, y))\$mesh")))
end

function r_ahull_arcs(points, alpha)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y alpha
    RCall.reval("library(alphahull)")
    return Float64.(RCall.rcopy(RCall.reval("alphahull::ahull(cbind(x, y), alpha=alpha)\$arcs")))
end

function prototype_edges(points)
    triangulation = RangeBuilder._shull_topology_prototype(points)
    return Set((min(edge[1], edge[2]), max(edge[1], edge[2]))
               for edge in RangeBuilder.DelaunayTriangulation.each_solid_edge(triangulation))
end

function final_shull_edges(points)
    state = RangeBuilder._shull_final_triads(points)
    return Set((min(a, b), max(a, b))
               for triad in state.triads
               for (a, b) in ((triad.a, triad.b), (triad.b, triad.c), (triad.a, triad.c)))
end

@testset "SHull Float32 topology prototype" begin
    rng = MersenneTwister(20260815)
    for point_count in (3, 4, 5, 7, 12, 24, 48, 96)
        points = rand(rng, Float64, point_count, 2) .* [320.0 160.0] .- [160.0 80.0]
        @test prototype_edges(points) == r_shull_edges(points)
        @test final_shull_edges(points) == r_shull_edges(points)
        state = RangeBuilder._shull_incremental_triads(points)
        @test Set(point.id for point in state.hull) == r_shull_hull(points)
        @test RangeBuilder._shull_trlist(points) == r_shull_trlist(points)
        @test isapprox(RangeBuilder._shull_mesh(points), r_delvor_mesh(points); rtol=2e-6, atol=1e-10)
    end

    alpha_points = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.4 0.5]
    alpha = 0.7
    julia_hull = RangeBuilder.ahull(RangeBuilder.delvor(alpha_points; backend=:shull); alpha)
    @test isapprox(julia_hull.arcs, r_ahull_arcs(alpha_points, alpha); rtol=2e-6, atol=1e-10)
end
