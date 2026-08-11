using RangeBuilder
using Random

function measure(label, f; repeats=3)
    f() # Compile and warm caches before measuring.
    GC.gc()
    samples = [@timed f() for _ in 1:repeats]
    order = sortperm(getfield.(samples, :time))
    sample = samples[order[cld(repeats, 2)]]
    println(rpad(label, 30), " ", round(sample.time; digits=4), " s   ",
            round(sample.bytes / 2.0^20; digits=2), " MiB allocated")
    return sample.value
end

n = 100_000
rng = MersenneTwister(20260809)
points = rand(rng, n, 2)
println("n = ", n, "; distribution = Uniform([0, 1]^2); seed = 20260809")

dv = measure("delvor", () -> delvor(points))
println(rpad("Delaunay/Voronoi mesh rows", 30), size(dv.mesh, 1))

for alpha in (0.01, 0.07)
    println("alpha = ", alpha)
    shape = measure("ashape (cached delvor)", () -> ashape(dv; alpha))
    hull = measure("ahull (cached delvor)", () -> ahull(dv; alpha))
    println(rpad("alpha-shape edges", 30), size(shape.edges, 1))
    println(rpad("alpha-hull arcs", 30), size(hull.arcs, 1))
    measure("areaahull (cached hull)", () -> areaahull(hull))
end

measure("ahull + areaahull alpha=0.07", () -> begin
    hull = ahull(points; alpha=0.07)
    areaahull(hull)
end)
