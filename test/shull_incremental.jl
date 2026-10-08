function _shull_edge_key(a::Int, b::Int)
    return a < b ? (a, b) : (b, a)
end

function _assert_shull_incremental_invariants(state, point_count::Int)
    triads = state.triads
    edge_adjacency = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for (index, triad) in enumerate(triads)
        @test length(unique((triad.a, triad.b, triad.c))) == 3
        for (edge, neighbour) in (
            (_shull_edge_key(triad.a, triad.b), triad.ab),
            (_shull_edge_key(triad.b, triad.c), triad.bc),
            (_shull_edge_key(triad.a, triad.c), triad.ac),
        )
            push!(get!(edge_adjacency, edge, Tuple{Int, Int}[]), (index, neighbour))
        end
    end

    @test length(triads) == 2point_count - length(state.hull) - 2
    for references in values(edge_adjacency)
        if length(references) == 1
            @test references[1][2] == 0
        else
            @test length(references) == 2
            first, second = references
            @test first[2] == second[1]
            @test second[2] == first[1]
        end
    end
end

@testset "SHull incremental topology" begin
    rng = MersenneTwister(20260815)
    for point_count in 3:24
        points = rand(rng, point_count, 2)
        state = RangeBuilder._shull_incremental_triads(points)
        @test length(state.points) == point_count
        @test length(unique(point.id for point in state.hull)) == length(state.hull)
        _assert_shull_incremental_invariants(state, point_count)

        final_state = RangeBuilder._shull_final_triads(points)
        _assert_shull_incremental_invariants(final_state, point_count)
    end
end
