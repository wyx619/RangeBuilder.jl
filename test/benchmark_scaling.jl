using RangeBuilder
using Random

function measure(label, f; repeats=3)
    f() # Compile and warm caches before measuring.
    GC.gc()
    samples = [@timed f() for _ in 1:repeats]
    order = sortperm(getfield.(samples, :time))
    sample = samples[order[cld(repeats, 2)]]
    println(rpad(label, 34), " ", round(sample.time; digits=4), " s   ",
            round(sample.bytes / 2.0^20; digits=2), " MiB allocated")
    return sample.value
end

for n in (4_000, 8_000, 10_000, 20_000)
    rng = MersenneTwister(20260809 + n)
    points = rand(rng, n, 2)
    println("\nn = ", n, "; alpha = 0.07; seed = ", 20260809 + n)
    dv = measure("delvor", () -> delvor(points))
    println(rpad("Delaunay/Voronoi mesh rows", 34), size(dv.mesh, 1))
    hull = measure("ahull (cached delvor)", () -> ahull(dv; alpha=0.07))
    println(rpad("alpha-hull arcs", 34), size(hull.arcs, 1))
    measure("areaahull (cached hull)", () -> areaahull(hull))
    measure("ahull + areaahull end-to-end", () -> begin
        h = ahull(points; alpha=0.07)
        areaahull(h)
    end)
end
