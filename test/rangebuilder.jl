@testset "RangeBuilder workflow" begin
    close_points = [0.0 0.0; 0.001 0.0; 1.0 1.0]
    @test filterByProximity(close_points, 1.0; returnIndex=true) == [1]
    @test size(filterByProximity(close_points, 1.0), 1) == 2
    @test ismissing(filterByProximity([0.0 0.0; 2.0 0.0], 1.0; returnIndex=true))

    synthetic_land = GeoInterface.Wrappers.Polygon([[(-1.0, -1.0), (1.0, -1.0),
                                                      (1.0, 1.0), (-1.0, 1.0),
                                                      (-1.0, -1.0)]]; crs=4326)
    RangeBuilder._range_world_cache[50] = synthetic_land
    land_points = [0.0 0.0; 10.0 0.0; NaN 0.0]
    land_mask = filterByLand(land_points; coastScale=50)
    @test land_mask[1] === true
    @test land_mask[2] === false
    @test ismissing(land_mask[3])

    empty!(RangeBuilder._range_world_cache)
    empty!(RangeBuilder._range_land_buffer_cache)
    @test filterByLand([116.0 39.0; 0.0 0.0]; coastScale=50) == Union{Missing,Bool}[true, false]

    antimeridian_points = [
        -179.7058 71.07439; -179.7584 71.1759; -179.4098 71.17809;
        -173.141667 64.78; -173.086667 64.83; 174.33 65.5;
        -179.583 70.95; -171.226 65.645
    ]
    antimeridian_range = getDynamicAlphaHull(
        antimeridian_points;
        fraction=0.95,
        partCount=3,
        buff=10_000,
        initialAlpha=2,
        alphaIncrement=1,
        alphaCap=400,
        clipToCoast=:terrestrial,
        backend=:shull,
        rng=RangeBuilder.RSeed.r_rng([
            10407, -1869255335, 1054321574, 982794370,
            -306250840, -2080751448, -1300630003,
        ]),
    )
    @test antimeridian_range.alpha == "alphaMCH"
    @test antimeridian_range.hull !== nothing

    dynamic_points = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.5]
    dynamic = getDynamicAlphaHull(dynamic_points; fraction=0.8, partCount=3, buff=0,
                                  initialAlpha=0.7, alphaIncrement=0.1,
                                  alphaCap=2, clipToCoast=:no)
    # The original R implementation still projects and runs GEOS buffer(0),
    # whose normalized boundary does not retain the required coverage here.
    @test dynamic.alpha == "alphaMCH"
    @test occursin("Polygon", string(typeof(dynamic.hull)))
    dynamic_shull = getDynamicAlphaHull(dynamic_points; fraction=0.8, partCount=3, buff=0,
                                        initialAlpha=0.7, alphaIncrement=0.1,
                                        alphaCap=2, clipToCoast=:no, backend=:shull)
    # The strict SHull path mirrors original R `interp` triad order, while
    # the default JuliaGeometry path intentionally retains its faster native
    # Delaunay traversal. They need not select the same alpha.
    @test dynamic_shull.alpha == "alphaMCH"
    @test occursin("Polygon", string(typeof(dynamic_shull.hull)))
    dynamic_shull_r_rng = getDynamicAlphaHull(
        dynamic_points;
        fraction=0.8,
        partCount=3,
        buff=0,
        initialAlpha=0.7,
        alphaIncrement=0.1,
        alphaCap=2,
        clipToCoast=:no,
        backend=:shull,
        rng=RangeBuilder.RSeed.RMersenneTwister(1),
    )
    @test dynamic_shull_r_rng.alpha == "alphaMCH"
    dynamic_buffered = getDynamicAlphaHull(dynamic_points; fraction=0.8, partCount=3, buff=1000,
                                           initialAlpha=0.7, alphaIncrement=0.1,
                                           alphaCap=2, clipToCoast=:no)
    @test dynamic_buffered.alpha == "alpha0.8"

    buffer_precision_source = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (0.01, 0.0),
                                                               (0.01, 0.01), (0.0, 0.01),
                                                               (0.0, 0.0)]]; crs=4326)
    buffer_precision = RangeBuilder._buffer_range(buffer_precision_source, 1000; force=true)
    @test GeoInterface.npoint(buffer_precision) == 125

    buffer_components = [
        GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (0.01, 0.0),
                                        (0.01, 0.01), (0.0, 0.01),
                                        (0.0, 0.0)]]; crs=4326),
        GeoInterface.Wrappers.Polygon([[(0.015, 0.0), (0.025, 0.0),
                                        (0.025, 0.01), (0.015, 0.01),
                                        (0.015, 0.0)]]; crs=4326),
    ]
    buffer_source = GeoInterface.Wrappers.MultiPolygon(buffer_components; crs=4326)
    buffered_components = RangeBuilder._buffer_range(buffer_source, 1000; force=true)
    @test GeoInterface.geomtrait(buffered_components) isa GeoInterface.GeometryCollectionTrait
    @test RangeBuilder._polygon_count(buffered_components) == 2
    @test startswith(RangeBuilder._geometry_wkt(buffered_components), "GEOMETRYCOLLECTION")

    dynamic_wkt = getDynamicAlphaHullWKT(dynamic_points; fraction=0.8, partCount=3, buff=0,
                                         initialAlpha=0.7, alphaIncrement=0.1,
                                         alphaCap=2, clipToCoast=:no)
    @test startswith(dynamic_wkt, "POLYGON") || startswith(dynamic_wkt, "MULTIPOLYGON")
    mch = getDynamicAlphaHull(dynamic_points; fraction=1, partCount=1, buff=0,
                              initialAlpha=0.1, alphaIncrement=0.1,
                              alphaCap=0.1, clipToCoast=:no)
    @test mch.alpha == "alphaMCH"

    batch_result = buildRanges(Dict(
        :valid => dynamic_points,
        :one_point => [0.0 0.0],
        :two_points => [0.0 0.0; 1.0 1.0],
        :collinear => [0.0 0.0; 1.0 1.0; 2.0 2.0],
    ); fraction=0.8, partCount=3, buff=0, initialAlpha=0.7,
       alphaIncrement=0.1, alphaCap=2, clipToCoast=:no)
    @test Set(keys(batch_result.ranges)) == Set((:valid,))
    @test batch_result.excluded[:one_point] == (status=:insufficient_points, npoints=1)
    @test batch_result.excluded[:two_points] == (status=:insufficient_points, npoints=2)
    @test batch_result.excluded[:collinear] == (status=:degenerate_geometry, npoints=3)

    square = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (1.0, 0.0),
                                              (1.0, 1.0), (0.0, 1.0),
                                              (0.0, 0.0)]]; crs=4326)
    stack = rasterStackFromPolyList((sp_a=square, sp_b=square); resolution=0.5)
    @test stack isa Rasters.RasterStack
    @test keys(stack) == (:sp_a, :sp_b)
    @test size(stack[:sp_a]) == (2, 2)
    @test all(parent(stack[:sp_a]) .== 1)
    @test Rasters.extent(stack[:sp_a]) == Rasters.Extent(X=(0.0, 1.0), Y=(0.0, 1.0))
    richness = speciesRichness(stack)
    @test size(richness) == (2, 2)
    @test all(parent(richness) .== 2)
    @test speciesRichness(stack; zeroToMissing=false)[1, 1] == 2

    expanded = rasterStackFromPolyList((sp_a=square, sp_b=square);
                                        resolution=1.0, extent=[0.0, 2.0, 0.0, 2.0])
    expanded_richness = speciesRichness(expanded)
    @test count(ismissing, parent(expanded_richness)) == 3
    @test count(==(2), skipmissing(parent(expanded_richness))) == 1

    tiny = GeoInterface.Wrappers.Polygon([[(0.10, 0.10), (0.20, 0.10),
                                            (0.20, 0.20), (0.10, 0.20),
                                            (0.10, 0.10)]]; crs=4326)
    tiny_stack = rasterStackFromPolyList(Dict("tiny" => tiny); resolution=1.0)
    @test parent(tiny_stack[:tiny])[1, 1] == 1
    mixed_stack = rasterStackFromPolyList((square=square, tiny=tiny);
                                           resolution=0.5, retainSmallRanges=false)
    @test keys(mixed_stack) == (:square,)
end
