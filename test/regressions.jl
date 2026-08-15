@testset "Regression cases" begin
    # Historical input preprocessing cases.
    collinear_rows = [0.0 0.0; 1.0 0.0; 2.0 0.0; 0.0 1.0]
    shuffled_rows = RangeBuilder._range_shuffle_collinear(collinear_rows; rng=MersenneTwister(0))
    @test !(shuffled_rows[1, 1] == shuffled_rows[2, 1] == shuffled_rows[3, 1])
    @test !(shuffled_rows[1, 2] == shuffled_rows[2, 2] == shuffled_rows[3, 2])
    @test Set(Tuple.(eachrow(shuffled_rows))) == Set(Tuple.(eachrow(collinear_rows)))
    collinear_cols = [0.0 0.0; 0.0 1.0; 0.0 2.0; 1.0 0.0]
    shuffled_cols = RangeBuilder._range_shuffle_collinear(collinear_cols; rng=MersenneTwister(0))
    @test !(shuffled_cols[1, 1] == shuffled_cols[2, 1] == shuffled_cols[3, 1])
    @test !(shuffled_cols[1, 2] == shuffled_cols[2, 2] == shuffled_cols[3, 2])
    @test RangeBuilder._range_shuffle_collinear([0.0 0.0; 1.0 1.0; 2.0 0.0]) ==
          [0.0 0.0; 1.0 1.0; 2.0 0.0]

    exact_dup_input = [0.0 0.0; 0.0 0.0; 1.0 1.0; 0.0 2.0]
    dropped_dv = RangeBuilder._range_drop_duplicate_points(exact_dup_input)
    @test dropped_dv !== nothing
    @test size(dropped_dv.x, 1) == 3
    @test RangeBuilder._range_drop_duplicate_points([0.0 0.0; 1e-10 1e-10; 2e-10 2e-10]) === nothing
    @test RangeBuilder._range_drop_duplicate_points([0.0 0.0; 1e-13 0.0; 1.0 1.0; 0.0 2.0]) !== nothing
    float32_collision = [116.0 39.0; 116.0000001 39.0; 117.0 39.0;
                         116.0 40.0; 117.0 40.0]
    @test RangeBuilder._shull_float32_duplicate_pair(float32_collision) == (1, 2)
    shull_dropped = RangeBuilder._range_drop_duplicate_points(float32_collision; backend=:shull)
    @test shull_dropped !== nothing
    @test size(shull_dropped.x, 1) == 4

    valid_points = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.4 0.5]
    valid_poly = ah2polygon(ahull(valid_points; alpha=0.7))
    @test RangeBuilder._range_polygon_is_valid(valid_poly)
    @test RangeBuilder._range_projected_polygon_is_valid(valid_poly)
    bowtie = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (1.0, 1.0), (1.0, 0.0),
                                             (0.0, 1.0), (0.0, 0.0)]]; crs=4326)
    @test !RangeBuilder._range_polygon_is_valid(bowtie)
    @test !RangeBuilder._range_projected_polygon_is_valid(bowtie)
    @test !RangeBuilder._range_polygon_is_valid(nothing)
    @test !RangeBuilder._range_projected_polygon_is_valid(nothing)

    holey_poly = GeoInterface.Wrappers.Polygon(
        [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
         [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]]; crs=4326)
    holey_pts = [3.0 5.0; 5.0 5.0; 5.0 0.0; 15.0 5.0; 0.0 5.0; 3.0 3.0]
    @test RangeBuilder._range_points_in_polygon(holey_pts, holey_poly) == [true, false, true, false, true, true]
    @test RangeBuilder._range_points_in_polygon(holey_pts, nothing) == falses(6)

    overlap_a = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (2.0, 0.0),
                                                 (2.0, 1.0), (0.0, 1.0),
                                                 (0.0, 0.0)]]; crs=4326)
    overlap_b = GeoInterface.Wrappers.Polygon([[(1.0, 0.0), (3.0, 0.0),
                                                 (3.0, 1.0), (1.0, 1.0),
                                                 (1.0, 0.0)]]; crs=4326)
    overlap = GeoInterface.Wrappers.GeometryCollection([overlap_a, overlap_b]; crs=4326)
    overlap_pts = [0.5 0.5; 1.5 0.5; 2.5 0.5]
    @test RangeBuilder._range_points_in_polygon(overlap_pts, overlap) == trues(3)
    @test RangeBuilder._range_point_intersection_counts(overlap_pts, overlap) == [1, 2, 1]

    proj_refs = [(0.0, 0.0, 0.0, 0.0),
                 (90.0, 45.0, 7396237.3744978, 5466867.76021372),
                 (-120.0, -30.0, -10758407.9263306, -3764325.42689213),
                 (179.0, 80.0, 10585770.1667077, 8205603.09648845),
                 (45.0, 0.0, 4310989.76555424, 0.0)]
    proj_ctx = RangeBuilder._range_native_context()
    for (lon, lat, rx, ry) in proj_refs
        px, py = proj_ctx.forward(lon, lat)
        @test isapprox(px, rx; rtol=1e-7)
        @test isapprox(py, ry; rtol=1e-7)
    end

    # Cocircular and near-collinear sets can have more than one valid Delaunay
    # triangulation.  Julia must be deterministic and return valid geometry;
    # their internal mesh and arc row order are intentionally not compared to R.
    cocircular = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.5]
    @test delvor(cocircular).mesh == delvor(cocircular).mesh
    cocircular_range = getDynamicAlphaHull(cocircular; fraction=0.8, partCount=3,
                                           buff=1000, initialAlpha=0.7,
                                           alphaIncrement=0.1, alphaCap=2,
                                           clipToCoast=:no)
    @test cocircular_range.alpha == "alpha0.8"
    @test RangeBuilder._range_polygon_is_valid(cocircular_range.hull)

    near_collinear = [0.0 0.0; 1.0 1e-9; 2.0 -1e-9; 3.0 2e-9; 1.5 0.4]
    near_hull = ahull(near_collinear; alpha=0.45)
    near_hull_again = ahull(near_collinear; alpha=0.45)
    @test size(near_hull.arcs, 2) == 8
    @test near_hull.xahull == near_hull_again.xahull
    @test near_hull.arcs == near_hull_again.arcs

    @test getDynamicAlphaHull([0.0 0.0; 1.0 1.0; 2.0 2.0]; clipToCoast=:no).alpha == "alphaMCH"
    @test_throws ArgumentError getDynamicAlphaHull([0.0 0.0; 0.0 1.0]; clipToCoast=:no)
end
