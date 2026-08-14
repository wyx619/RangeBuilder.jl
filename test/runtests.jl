using Test
using RangeBuilder
using Random
using Rasters
import GeoInterface

pts = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.4 0.5]
dv = delvor(pts)
@test size(dv.x) == (5, 2)
@test size(dv.mesh, 2) == 12
@test size(dv.mesh, 1) >= 5
@test all(isfinite, dv.mesh[:, 3:10])
@test RangeBuilder.DelaunayTriangulation.has_ghost_vertices(
    RangeBuilder.DelaunayTriangulation.get_graph(dv.triangulation)
)

sh = ashape(dv; alpha=0.5)
@test sh.delvor === dv
@test sh.alpha == 0.5
@test size(sh.edges, 2) == 12
@test sh.length >= 0

cp = complement(dv; alpha=0.5)
@test size(cp, 2) == 19
@test all(isfinite, cp[:, 1:2])

ah = ahull(dv; alpha=0.5)
@test size(ah.arcs, 2) == 8
@test ah.length >= 0
@test ah.alpha == 0.5
poly = ah2polygon(ahull(pts; alpha=0.7))
@test poly !== nothing
@test occursin("Polygon", string(typeof(poly)))
@test inahull(ah, [2.0 2.0]) == [false]
@test eltype(inahull(ah, [0.5 0.5])) == Bool

ci = inter(0, 0, 1, 1, 0, 1)
@test ci.n_cut == 2
@test ci.v1 == (1.0, 0.0)
@test ci.theta1 ≈ pi / 3
@test inter(0, 0, 1, 3, 0, 1).n_cut == 0
@test areaahull(ah) >= 0
@test rotation([1.0, 0.0], pi / 2) ≈ [0.0, -1.0]
@test anglesArc([1.0, 0.0], pi / 4) ≈ [-pi / 4, pi / 4]
@test size(koch(side=1, niter=3), 1) == 48
@test size(rkoch(10; side=1, niter=2, rng=MersenneTwister(1)), 1) == 10

wpts = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
w = dw(wpts; eps=0.7)
@test size(w, 2) == 9
@test all(w[:, 3] .== 0.7)
@test size(arc(w[1, 1:2], w[1, 3], w[1, 4:5], w[1, 6])) == (100, 2)
@test hasmethod(dw_track, Tuple{AbstractMatrix})
@test hasmethod(ahull_track, Tuple{AbstractMatrix})

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
bundled_land_mask = filterByLand([116.0 39.0; 0.0 0.0]; coastScale=50)
@test bundled_land_mask == Union{Missing,Bool}[true, false]

dynamic = getDynamicAlphaHull(pts; fraction=0.8, partCount=3, buff=0,
                              initialAlpha=0.7, alphaIncrement=0.1,
                              alphaCap=2, clipToCoast=:no)
@test dynamic.alpha == "alpha0.8"
@test occursin("Polygon", string(typeof(dynamic.hull)))
dynamic_buffered = getDynamicAlphaHull(pts; fraction=0.8, partCount=3, buff=1000,
                                       initialAlpha=0.7, alphaIncrement=0.1,
                                       alphaCap=2, clipToCoast=:no)
@test dynamic_buffered.alpha == "alpha0.8"

buffer_precision_source = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (0.01, 0.0),
                                                           (0.01, 0.01), (0.0, 0.01),
                                                           (0.0, 0.0)]]; crs=4326)
buffer_precision = RangeBuilder._buffer_range(buffer_precision_source, 1000; force=true)
@test GeoInterface.npoint(buffer_precision) == 45

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

dynamic_wkt = getDynamicAlphaHullWKT(pts; fraction=0.8, partCount=3, buff=0,
                                     initialAlpha=0.7, alphaIncrement=0.1,
                                     alphaCap=2, clipToCoast=:no)
@test startswith(dynamic_wkt, "POLYGON") || startswith(dynamic_wkt, "MULTIPOLYGON")
mch = getDynamicAlphaHull(pts; fraction=1, partCount=1, buff=0,
                          initialAlpha=0.1, alphaIncrement=0.1,
                          alphaCap=0.1, clipToCoast=:no)
@test mch.alpha == "alphaMCH"

batch_result = buildRanges(Dict(
    :valid => pts,
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

@test_throws ArgumentError delvor([0.0 0.0; 1.0 1.0; 2.0 2.0])
@test_throws ArgumentError delvor([0.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 0.0])
@test_throws ArgumentError ashape(pts; alpha=-1)
@test_throws ArgumentError dw(pts; eps=0)

@test occursin("DelVor", sprint(show, dv))
@test occursin("mesh", sprint(show, MIME"text/plain"(), dv))
@test occursin("AShape", sprint(show, sh))
@test occursin("alpha_extremes", sprint(show, MIME"text/plain"(), sh))
@test occursin("AHull", sprint(show, ah))
@test occursin("complement", sprint(show, MIME"text/plain"(), ah))
@test occursin("CircleIntersection", sprint(show, ci))

# _range_shuffle_collinear: R-style collinear preprocessing (plan Task 1)
collinear_rows = [0.0 0.0; 1.0 0.0; 2.0 0.0; 0.0 1.0]  # first three rows share y == 0
shuffled_rows = RangeBuilder._range_shuffle_collinear(collinear_rows; rng=MersenneTwister(0))
@test !(shuffled_rows[1, 1] == shuffled_rows[2, 1] == shuffled_rows[3, 1])
@test !(shuffled_rows[1, 2] == shuffled_rows[2, 2] == shuffled_rows[3, 2])
@test Set(Tuple.(eachrow(shuffled_rows))) == Set(Tuple.(eachrow(collinear_rows)))
collinear_cols = [0.0 0.0; 0.0 1.0; 0.0 2.0; 1.0 0.0]  # first three rows share x == 0
shuffled_cols = RangeBuilder._range_shuffle_collinear(collinear_cols; rng=MersenneTwister(0))
@test !(shuffled_cols[1, 1] == shuffled_cols[2, 1] == shuffled_cols[3, 1])
@test !(shuffled_cols[1, 2] == shuffled_cols[2, 2] == shuffled_cols[3, 2])
non_collinear_rows = [0.0 0.0; 1.0 1.0; 2.0 0.0]
@test RangeBuilder._range_shuffle_collinear(non_collinear_rows) == non_collinear_rows
@test RangeBuilder._range_shuffle_collinear(pts; rng=MersenneTwister(0)) == pts

# _range_drop_duplicate_points: R-style duplicate-point drop (plan Task 2)
exact_dup_input = [0.0 0.0; 0.0 0.0; 1.0 1.0; 0.0 2.0]
dropped_dv = RangeBuilder._range_drop_duplicate_points(exact_dup_input)
@test dropped_dv !== nothing
@test size(dropped_dv.x, 1) == 3
@test RangeBuilder._range_drop_duplicate_points([0.0 0.0; 1e-10 1e-10; 2e-10 2e-10]) === nothing
near_noncollinear = [0.0 0.0; 1e-13 0.0; 1.0 1.0; 0.0 2.0]
@test RangeBuilder._range_drop_duplicate_points(near_noncollinear) !== nothing

# _range_polygon_is_valid: GEOS topological validity (plan Task 3)
valid_poly = RangeBuilder.ah2polygon(RangeBuilder.ahull(pts; alpha=0.7))
@test RangeBuilder._range_polygon_is_valid(valid_poly)
bowtie = GeoInterface.Wrappers.Polygon([[(0.0, 0.0), (1.0, 1.0), (1.0, 0.0),
                                         (0.0, 1.0), (0.0, 0.0)]]; crs=4326)
@test !RangeBuilder._range_polygon_is_valid(bowtie)
@test !RangeBuilder._range_polygon_is_valid(nothing)
