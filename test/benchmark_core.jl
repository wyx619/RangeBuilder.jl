using RangeBuilder
using Random

function measure(label, f; repeats=7)
    f() # Compile and warm caches before measuring.
    GC.gc()
    results = [@timed f() for _ in 1:repeats]
    times = sort([result.time for result in results])
    median_time = times[cld(length(times), 2)]
    median_result = results[sortperm([result.time for result in results])[cld(length(results), 2)]]
    println("$(rpad(label, 24))  $(round(median_time; digits=4)) s  " *
            "$(round(median_result.bytes / 2.0^20; digits=2)) MiB allocated")
    return median_result.value
end

for n in (100, 500, 2000)
    println("\nn = $n, alpha = 0.07")
    rng = MersenneTwister(20260808 + n)
    points = rand(rng, n, 2)

    dv = measure("delvor", () -> delvor(points))
    println("$(rpad("mesh edges", 24))  $(size(dv.mesh, 1))")
    measure("ashape (cached delvor)", () -> ashape(dv; alpha=0.07))
    measure("complement (cached)", () -> complement(dv; alpha=0.07))
    hull = measure("ahull (cached delvor)", () -> ahull(dv; alpha=0.07))
    println("$(rpad("alpha-hull arcs", 24))  $(size(hull.arcs, 1))")
    n == 100 && measure("dw", () -> dw(points; eps=0.07))
end
