using Random
using RCall

points = rand(MersenneTwister(20260811 + 10_000), 10_000, 2)
@rput points

RCall.reval("""
    library(rangeBuilder)
    .rangebuilder_once <- function(points) {
        getDynamicAlphaHull(
            points,
            fraction = 0.95,
            partCount = 3,
            buff = 0,
            coordHeaders = c(1, 2),
            clipToCoast = "no",
            initialAlpha = 0.005,
            alphaIncrement = 0.005,
            alphaCap = 0.050
        )
    }
    .rangebuilder_once_result <- NULL
    .rangebuilder_once_time <- unname(system.time(
        .rangebuilder_once_result <- .rangebuilder_once(points)
    )["elapsed"])
""")

elapsed = RCall.rcopy(RCall.reval(".rangebuilder_once_time"))
alpha = RCall.rcopy(RCall.reval(".rangebuilder_once_result[[2]]"))
println("R rangeBuilder dynamic search n=10000: $(elapsed) s; selected alpha=$(alpha)")
