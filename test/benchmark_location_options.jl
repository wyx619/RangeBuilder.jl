using DelaunayTriangulation
using Random

function measure(label, f; repeats=2)
    f()
    GC.gc()
    times = [@elapsed f() for _ in 1:repeats]
    println(rpad(label, 38), " ", round(sort(times)[cld(repeats, 2)]; digits=4), " s")
end

for n in (10_000, 20_000)
    points = rand(MersenneTwister(20260809 + n), n, 2)
    pts = [(points[i, 1], points[i, 2]) for i in axes(points, 1)]
    println("\nn = ", n)
    measure("default", () -> triangulate(pts; rng=MersenneTwister(7)))
    measure("num_sample_rule = _ -> 1", () -> triangulate(pts; rng=MersenneTwister(7), num_sample_rule = _ -> 1))
    measure("recompute_representative_points = false", () -> triangulate(pts; rng=MersenneTwister(7), recompute_representative_points=false))
    measure("check_arguments = false", () -> triangulate(pts; rng=MersenneTwister(7), check_arguments=false))
    measure("FastKernel", () -> triangulate(pts; rng=MersenneTwister(7), predicates=FastKernel()))
end
