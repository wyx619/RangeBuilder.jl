@testset "R-compatible RNG" begin
    rng = RangeBuilder.RSeed.RMersenneTwister(1)
    expected = [
        0.26550866314209998, 0.37212389963679016, 0.57285336335189641,
        0.90820778999477625, 0.2016819310374558, 0.89838968496769667,
        0.94467526860535145, 0.66079779248684645, 0.62911404389888048,
        0.061786270467564464,
    ]
    actual = [RangeBuilder.RSeed.r_unif(rng) for _ in eachindex(expected)]
    @test actual == expected

    rng = RangeBuilder.RSeed.RMersenneTwister(1)
    @test RangeBuilder.RSeed.r_runif(rng, 4, -1, 1) ==
          2 .* expected[1:4] .- 1

    rng = RangeBuilder.RSeed.RMersenneTwister(1)
    @test RangeBuilder.RSeed.r_jitter(rng, [10.0, 20.0], 0.5) ==
          [10.0 + 0.5 * (2 * expected[1] - 1),
           20.0 + 0.5 * (2 * expected[2] - 1)]

    source = RangeBuilder.RSeed.RMersenneTwister(123)
    RangeBuilder.RSeed.r_runif(source, 17)
    r_state = Vector{Int32}(undef, 626)
    r_state[1] = 10403
    r_state[2] = source.index
    for index in eachindex(source.state)
        r_state[index + 2] = reinterpret(Int32, source.state[index])
    end
    restored = RangeBuilder.RSeed.RMersenneTwister(r_state)
    @test RangeBuilder.RSeed.r_runif(restored, 10) == RangeBuilder.RSeed.r_runif(source, 10)

    invalid_kind = copy(r_state)
    invalid_kind[1] = 10402
    @test_throws ArgumentError RangeBuilder.RSeed.RMersenneTwister(invalid_kind)
    @test_throws ArgumentError RangeBuilder.RSeed.RMersenneTwister(r_state[1:end-1])

    lecuyer_expected = [
        0.6775328286287442, 0.42734572288764422, 0.9103805304875483,
        0.95572819835307676, 0.84065858527482162, 0.34366115612944603,
        0.35657673239893289, 0.42839009922583143, 0.44749898861157472,
        0.40259352669572779,
    ]
    lecuyer = RangeBuilder.RSeed.RLecuyerCMRG(1)
    @test [RangeBuilder.RSeed.r_unif(lecuyer) for _ in eachindex(lecuyer_expected)] == lecuyer_expected

    lecuyer_state = Int32[
        10407, 1280795612, -169270483, -442010614, -603558397,
        -222347416, 1489374793,
    ]
    restored_lecuyer = RangeBuilder.RSeed.r_rng(lecuyer_state)
    @test restored_lecuyer isa RangeBuilder.RSeed.RLecuyerCMRG
    @test RangeBuilder.RSeed.r_runif(restored_lecuyer, 10) == lecuyer_expected
    @test_throws ArgumentError RangeBuilder.RSeed.r_rng(Int32[10402, 1])

    points = [0.0 0.0; 1.0 10.0; 2.0 20.0]
    rng = RangeBuilder.RSeed.RMersenneTwister(1)
    jittered = RangeBuilder._shull_jittered_points(points; rng)
    @test jittered[:, 1] ==
          RangeBuilder.RSeed.r_jitter(RangeBuilder.RSeed.RMersenneTwister(1), points[:, 1], 0.0002)
    rng = RangeBuilder.RSeed.RMersenneTwister(1)
    RangeBuilder.RSeed.r_runif(rng, size(points, 1), -1.0, 1.0)
    @test jittered[:, 2] == RangeBuilder.RSeed.r_jitter(rng, points[:, 2], 0.002)
end
