# Optional development-only verification against R's live RNG implementation.
#
# Run explicitly with: julia --project=test test/rng_reference.jl
# RCall is a test-only dependency; RangeBuilder.jl remains pure Julia at runtime.
using RCall
using Test
using RangeBuilder

RCall.reval("options(warn = -1)")

@testset "Original R reference: RNG and jitter" begin
    RCall.reval(
        "RNGkind(\"L'Ecuyer-CMRG\"); " *
        "set.seed(20260816L); " *
        ".reference_rng_seed <- .Random.seed; " *
        ".reference_runif <- runif(12); " *
        ".reference_jitter_x <- jitter(c(-4.5, 0.0, 7.25), factor=0.001); " *
        ".reference_jitter_y <- jitter(c(10.0, 12.0, 14.0), factor=0.001)",
    )

    r_seed = Int32.(RCall.rcopy(RCall.reval(".reference_rng_seed")))
    r_runif_values = Float64.(RCall.rcopy(RCall.reval(".reference_runif")))
    r_jitter_x = Float64.(RCall.rcopy(RCall.reval(".reference_jitter_x")))
    r_jitter_y = Float64.(RCall.rcopy(RCall.reval(".reference_jitter_y")))

    rng = RangeBuilder.RSeed.r_rng(r_seed)
    @test RangeBuilder.RSeed.r_runif(rng, length(r_runif_values)) == r_runif_values

    # base::jitter() derives its default amount from the smallest rounded
    # adjacent spacing; `_shull_r_jitter_amount` reproduces that calculation.
    x = [-4.5, 0.0, 7.25]
    y = [10.0, 12.0, 14.0]
    @test RangeBuilder.RSeed.r_jitter(rng, x, RangeBuilder._shull_r_jitter_amount(x)) ==
          r_jitter_x
    @test RangeBuilder.RSeed.r_jitter(rng, y, RangeBuilder._shull_r_jitter_amount(y)) ==
          r_jitter_y
end
