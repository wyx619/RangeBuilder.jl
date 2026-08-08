using AlphaHull

function halton(i, base)
    value = 0.0
    factor = 1.0 / base
    while i > 0
        value += (i % base) * factor
        i = div(i, base)
        factor /= base
    end
    return value
end

function benchmark_points(n)
    points = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        points[i, 1] = halton(i, 2)
        points[i, 2] = halton(i, 3)
    end
    return points
end

function median_time(f; repeats=7)
    f()
    GC.gc()
    times = sort([@elapsed f() for _ in 1:repeats])
    return times[cld(length(times), 2)]
end

for n in (500, 2000)
    points = benchmark_points(n)
    alpha = 0.07
    elapsed = median_time(() -> begin
        hull = ahull(points; alpha)
        areaahull(hull)
    end)
    println("Julia ahull + areaahull n=$n: $(round(elapsed; digits=4)) s")
end
