using DelaunayTriangulation
using Random

function measure(label, f; repeats=2)
    f()
    GC.gc()
    times = [@elapsed f() for _ in 1:repeats]
    println(rpad(label, 34), " ", round(sort(times)[cld(repeats, 2)]; digits=4), " s")
end

function run_case(label, n)
    points = rand(MersenneTwister(20260809 + n), n, 2)
    pts = [(points[i, 1], points[i, 2]) for i in axes(points, 1)]
    measure(label, () -> triangulate(pts; rng=MersenneTwister(7)))
end

for n in (4_000, 8_000, 10_000, 20_000)
    run_case("original n=$n", n)
end

@eval DelaunayTriangulation begin
    function has_ghost_vertices(G::Graph{I}) where {I}
        return I(𝒢) in get_vertices(G)
    end
end

for n in (4_000, 8_000, 10_000, 20_000)
    run_case("sentinel fix n=$n", n)
end
