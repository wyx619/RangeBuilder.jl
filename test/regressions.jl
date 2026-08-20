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
    # R chooses the first row returned by `which(..., arr.ind=TRUE)` on a
    # symmetric distance matrix, which is the second point of this pair.
    @test RangeBuilder._closest_point_pair([0.0 0.0; 1.0 0.0; 3.0 0.0]) == (2, 1)
    @test RangeBuilder._range_geodesic_closest_point_pair(
        [0.0 0.0; 1.0 0.0; 3.0 0.0],
    ) == (2, 1)
    # Dot-product comparisons saturate for sub-metre angular separations and
    # can retain an earlier, farther pair. Haversine comparison must select
    # the later, genuinely closer pair, following R's column-major tie rule.
    submetre_pairs = [0.0 0.0;
                       1.1e-6 0.0;
                       114.963105 29.326488;
                       114.963106 29.326488]
    @test RangeBuilder._range_geodesic_closest_point_pair(submetre_pairs) == (4, 3)

    exact_dup_input = [0.0 0.0; 0.0 0.0; 1.0 1.0; 0.0 2.0]
    dropped_dv = RangeBuilder._range_drop_duplicate_points(exact_dup_input)
    @test dropped_dv !== nothing
    @test size(dropped_dv.x, 1) == 3
    @test RangeBuilder._range_drop_duplicate_points([0.0 0.0; 1e-10 1e-10; 2e-10 2e-10]) === nothing
    @test RangeBuilder._range_drop_duplicate_points([0.0 0.0; 1e-13 0.0; 1.0 1.0; 0.0 2.0]) !== nothing
    float32_collision = [116.0 39.0; 116.0000001 39.0; 117.0 39.0;
                         116.0 40.0; 117.0 40.0]
    @test RangeBuilder._shull_float32_duplicate_pair(float32_collision) == (1, 2)
    signed_zero_collision = [-0.0 0.0; 0.0 -0.0; 1.0 1.0]
    @test RangeBuilder._shull_float32_duplicate_pair(signed_zero_collision) == (1, 2)
    @test_throws ArgumentError delvor(signed_zero_collision; backend=:shull)
    shull_duplicate_status = RangeBuilder._range_try_ahull_with_duplicate_status(
        float32_collision,
        1.0,
        :shull,
        MersenneTwister(0),
    )
    @test isnothing(shull_duplicate_status.hull)
    @test shull_duplicate_status.duplicate
    shull_dropped = RangeBuilder._range_drop_duplicate_points(float32_collision; backend=:shull)
    @test shull_dropped !== nothing
    @test size(shull_dropped.x, 1) == 4
    shull_prepared = RangeBuilder._range_drop_duplicate_points_with_coordinates(
        float32_collision; backend=:shull)
    @test shull_prepared !== nothing
    @test size(shull_prepared.points, 1) == 4
    @test shull_prepared.points == shull_prepared.delaunay.x

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

    # S2 validates geographic rings along great-circle arcs.  This ring is
    # planar-valid but self-crosses on the sphere, so it must not terminate a
    # dynamic alpha search before the R/S2-valid candidate is reached.
    spherical_crossing = ah2polygon(ahull(
        delvor([33.67 41.4905; 46.7739 39.6442; 46.3703 39.6719]; backend=:shull);
        alpha=97,
    ))
    @test !RangeBuilder._range_polygon_is_valid(spherical_crossing)

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

    disjoint_b = GeoInterface.Wrappers.Polygon([[(4.0, 0.0), (5.0, 0.0),
                                                   (5.0, 1.0), (4.0, 1.0),
                                                   (4.0, 0.0)]]; crs=4326)
    disjoint = GeoInterface.Wrappers.GeometryCollection([overlap_a, disjoint_b]; crs=4326)
    disjoint_pts = [0.5 0.5; 4.5 0.5; 8.0 0.5]
    @test RangeBuilder._range_point_intersection_counts(disjoint_pts, disjoint) == [1, 1, 0]

    # rangeBuilder uses sf with S2 enabled for WGS84 point coverage.  The
    # great-circle edge between (120, 0) and (0, 60) contains this point on
    # the sphere, while the longitude/latitude planar triangle does not.
    spherical_poly = GeoInterface.Wrappers.Polygon(
        [[(0.0, 0.0), (120.0, 0.0), (0.0, 60.0), (0.0, 0.0)]]; crs=4326)
    @test RangeBuilder._range_point_intersection_counts([60.0 45.0], spherical_poly) == [1]

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
@testset "SHull NaN circumcentre propagation" begin
    points = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0; 0.4 0.5]
    reference = delvor(points; backend=:shull)
    nonfinite_mesh = copy(reference.mesh)
    nonfinite_mesh[1, 7] = NaN
    nonfinite = DelVor(reference.x, nonfinite_mesh, reference.triangulation)

    # R's complement() aborts after NA reaches `sum(vert)`.  This must be a
    # candidate failure, rather than a silently altered Julia alpha hull.
    @test_throws ArgumentError complement(nonfinite; alpha=0.5)
    @test RangeBuilder._range_try_ahull(nonfinite, 0.5) === nothing

    # `tri.mesh()` restores original coordinates after a jitter retry, but C++
    # has already exported the triad using the Float64 jitter input. A
    # Float32 SHull state alone cannot retain that orientation.
    collinear = [0.0 0.0; 1.0 0.0; 2.0 0.0]
    triad = RangeBuilder._SHullTriad(1, 2, 3)
    jittered_points = [
        RangeBuilder._SHullPoint(1, 0f0, 0f0),
        RangeBuilder._SHullPoint(2, 1f0, 0f0),
        RangeBuilder._SHullPoint(3, 2f0, 1f-3),
    ]
    jittered_xy = [0.0 0.0; 1.0 0.0; 2.0 1e-3]
    jittered_state = (; points=jittered_points, triads=[triad], hull=copy(jittered_points),
                      orientation_xy=jittered_xy)
    @test RangeBuilder._shull_trlist(collinear, jittered_state)[1, 1:3] == [1, 2, 3]

    flipped_xy = [0.0 0.0; 1.0 0.0; 2.0 -1e-3]
    flipped_state = (; points=jittered_points, triads=[triad], hull=copy(jittered_points),
                     orientation_xy=flipped_xy)
    @test RangeBuilder._shull_trlist(collinear, flipped_state)[1, 1:3] == [1, 3, 2]

    # R's unique(sort(round(...))) merges signed zero. Julia's unique! uses
    # isequal, so this must be normalised before deriving jitter amplitude.
    signed_zero_values = [-120.0, -0.04, 0.04, 0.1, 124.5]
    @test RangeBuilder._shull_r_jitter_amount(signed_zero_values) == 2e-5
end

@testset "Spherical MCH convex hull" begin
    points = [0.0 0.0; 90.0 0.0; 0.0 60.0; 20.0 15.0]
    indices = RangeBuilder._spherical_convex_hull_indices(points)
    @test Set(indices) == Set((1, 2, 3))
    polygon = RangeBuilder._spherical_convex_hull_polygon(points)
    ring = GeoInterface.getpoint(GeoInterface.getexterior(polygon))
    @test length(ring) == 4
    @test ring[1] == ring[end]
end

@testset "Bundled one-degree map resource" begin
    payload = RangeBuilder._read_one_degree_map_resource()
    @test isfile(RangeBuilder._one_degree_map_path)
    @test startswith(payload["crs"], "+proj=cea")
    @test length(payload["ids"]) == length(payload["wkb"])
    @test length(payload["ids"]) > 10_000
    @test all(id -> id >= 0, payload["ids"])
    @test length(unique(payload["ids"])) == length(payload["ids"])
    @test all(bytes -> !isempty(bytes), payload["wkb"])
    first_geometry = RangeBuilder.LibGEOS.readgeom(
        Vector{Cuchar}(payload["wkb"][1]),
    )
    last_geometry = RangeBuilder.LibGEOS.readgeom(
        Vector{Cuchar}(payload["wkb"][end]),
    )
    @test first_geometry !== nothing
    @test last_geometry !== nothing
end
