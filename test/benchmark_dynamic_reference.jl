using Random
using RCall
using RangeBuilder

function benchmark_points(n)
    rng = MersenneTwister(20260811 + n)
    return rand(rng, n, 2)
end

function median_time(f; repeats=7)
    f()
    GC.gc()
    times = sort([@elapsed f() for _ in 1:repeats])
    return times[cld(length(times), 2)]
end

function r_dynamic_time(points; initialAlpha, alphaIncrement, alphaCap)
    @rput points
    RCall.reval("""
        library(rangeBuilder)
        .rangebuilder_dynamic <- function(points) {
            getDynamicAlphaHull(
                points,
                fraction = 0.95,
                partCount = 3,
                buff = 0,
                coordHeaders = c(1, 2),
                clipToCoast = "no",
                initialAlpha = $initialAlpha,
                alphaIncrement = $alphaIncrement,
                alphaCap = $alphaCap
            )
        }
        .rangebuilder_dynamic(points)
        .rangebuilder_times <- vapply(
            seq_len(7),
            function(i) unname(system.time(.rangebuilder_dynamic(points))["elapsed"]),
            numeric(1)
        )
        .rangebuilder_result <- .rangebuilder_dynamic(points)
    """)
    return (
        time=RCall.rcopy(RCall.reval("median(.rangebuilder_times)")),
        alpha=RCall.rcopy(RCall.reval(".rangebuilder_result[[2]]")),
    )
end

for (n, initialAlpha, alphaIncrement, alphaCap) in (
    (500, 0.005, 0.005, 0.020),
    (1_000, 0.005, 0.005, 0.050),
    (10_000, 0.005, 0.005, 0.050),
)
    points = benchmark_points(n)
    julia_result = getDynamicAlphaHull(
        points;
        fraction=0.95,
        partCount=3,
        buff=0,
        initialAlpha,
        alphaIncrement,
        alphaCap,
        clipToCoast=:no,
    )
    julia_time = median_time(() -> getDynamicAlphaHull(
        points;
        fraction=0.95,
        partCount=3,
        buff=0,
        initialAlpha,
        alphaIncrement,
        alphaCap,
        clipToCoast=:no,
    ))
    r_result = r_dynamic_time(points; initialAlpha, alphaIncrement, alphaCap)
    candidates = Int(round((alphaCap - initialAlpha) / alphaIncrement)) + 1
    println("dynamic search n=$n; candidates=$candidates; selected alpha: Julia=$(julia_result.alpha), R=$(r_result.alpha)")
    println("Julia getDynamicAlphaHull n=$n: $(round(julia_time; digits=4)) s")
    println("R rangeBuilder getDynamicAlphaHull n=$n: $(round(r_result.time; digits=4)) s")
end
