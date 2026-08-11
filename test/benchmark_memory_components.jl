using RangeBuilder
using DelaunayTriangulation
using Random

function measure(label, f; repeats=2)
    f()
    GC.gc()
    samples = [@timed f() for _ in 1:repeats]
    order = sortperm(getfield.(samples, :time))
    sample = samples[order[cld(repeats, 2)]]
    println(rpad(label, 34), " ", round(sample.time; digits=4), " s   ",
            round(sample.bytes / 2.0^20; digits=2), " MiB allocated")
    return sample.value
end

for n in (20_000, 100_000)
    points = rand(MersenneTwister(20260809 + n), n, 2)
    pts = [(points[i, 1], points[i, 2]) for i in axes(points, 1)]
    println("\nn = ", n)
    xy = measure("_points / Matrix conversion", () -> RangeBuilder._points(points))
    measure("point tuple conversion", () -> [(points[i, 1], points[i, 2]) for i in axes(points, 1)])
    tri = measure("triangulate only", () -> triangulate(pts; rng=MersenneTwister(7)))
    dv = measure("delvor full", () -> delvor(points))
    println(rpad("triangulation resident size", 34), " ",
            round(Base.summarysize(dv.triangulation) / 2.0^20; digits=2), " MiB")
    println(rpad("DelVor resident size", 34), " ",
            round(Base.summarysize(dv) / 2.0^20; digits=2), " MiB")
    hull = measure("ahull cached", () -> ahull(dv; alpha=0.07))
    measure("areaahull cached", () -> areaahull(hull))
end
