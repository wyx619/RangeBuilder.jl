@testset "Bundled optimized workflow" begin
    table = (
        Accepted_species=[:oak, :oak, :pine],
        Accepted_family=[:family_a, :family_a, :family_b],
        mylongitude=[116.0, 116.5, 0.0],
        mylatitude=[39.0, 39.4, 0.0],
    )
    input = standardize_occurrences(
        table;
        species=:Accepted_species,
        longitude=:mylongitude,
        latitude=:mylatitude,
        family=:Accepted_family,
        order=nothing,
    )
    @test input isa OccurrenceInput
    @test input.species == ["oak", "oak", "pine"]
    @test input.metadata.family == [:family_a, :family_a, :family_b]
    @test all(ismissing, input.metadata.order)
    @test input.source_rows == [1, 2, 3]

    rounded = standardize_occurrences(
        (species=["oak", "oak", "oak"],
         longitude=[116.011, 116.014, 116.061],
         latitude=[39.001, 39.004, 39.061]);
        coordinate_digits=2,
    )
    @test rounded.source_rows == [1, 3]
    @test rounded.longitude == [116.01, 116.06]
    @test rounded.latitude == [39.0, 39.06]
    one_digit = standardize_occurrences(
        (species=["oak", "oak"], longitude=[116.01, 116.04], latitude=[39.01, 39.04]);
        coordinate_digits=1,
    )
    @test one_digit.source_rows == [1]
    @test one_digit.longitude == [116.0]
    @test one_digit.latitude == [39.0]
    @test_throws ArgumentError standardize_occurrences(
        (species=["oak"], longitude=[1.0], latitude=[1.0]);
        coordinate_digits=-1,
    )

    grid = read_one_degree_map()
    @test length(grid.ids) == length(grid.geometries)
    @test length(grid.ids) > 10_000

    context = new_range_context(grid; buffer_m=10_000, clip_to_coast=:terrestrial)
    @test_throws ArgumentError read_one_degree_map(target_crs="EPSG:4326")
    @test_throws ErrorException new_range_context(grid; candidate_cache_limit=-1)
    points = Dict(
        :inland => [116.0 39.0; 116.5 39.0; 116.4 39.4; 116.1 39.3],
        :direct => [0.0 0.0; 0.1 0.1],
    )
    result = build_species_ranges(points, context; workers=1)

    @test length(result.results) == 2
    @test all(item -> item.success, result.results)
    @test all(item -> !isempty(item.grid_ids), result.results)
    @test result.results[1].species == :direct
    @test result.results[2].species == :inland
    @test result.s2_stats[:interior_covers] >= 1
    @test result.elapsed_seconds >= 0

    table_result = build_species_ranges(input, context; workers=1)
    oak = only(filter(item -> item.species == :oak, table_result.results))
    @test oak.metadata.family == :family_a
    @test ismissing(oak.metadata.order)
    @test oak.source_rows == [1, 2]

    degenerate = build_species_ranges(
        Dict(:line => [0.0 0.0; 1.0 1.0; 2.0 2.0]), context; workers=1)
    @test degenerate.results[1].method ==
          Symbol("alpha_hull_10.0km_buffer_clipped_to_land")
    @test degenerate.results[1].status != :degenerate_geometry

    float32_duplicate = build_species_ranges(
        Dict(:near_duplicate => [10.0 20.0; 10.0 + 1e-8 20.0 + 1e-8; 11.0 20.5]),
        context;
        workers=1,
    )
    @test float32_duplicate.results[1].success
    @test float32_duplicate.results[1].method == Symbol(
        "direct_point_grid_plus_10.0km_buffer_area_shull_precision_fallback")

    bundled_context = prepare_workflow_context(candidate_cache_limit=0)
    @test bundled_context.candidate_cache_limit == 0
end
