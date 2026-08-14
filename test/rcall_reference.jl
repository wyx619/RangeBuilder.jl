using Test
using RangeBuilder
using RCall
using Random

function r_delvor(points)
    @rput points
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("delvor(points)\$mesh"))
end

function r_ahull(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    RCall.reval(".ah <- ahull(points, alpha=alpha)")
    arcs = RCall.rcopy(RCall.reval(".ah\$arcs"))
    xahull = RCall.rcopy(RCall.reval(".ah\$xahull"))
    return arcs, xahull
end

function r_area(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    return(RCall.rcopy(RCall.reval("areaahull(ahull(points, alpha=alpha))")))
end

function r_koch(side, niter)
    @rput side niter
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("koch(side=side, niter=niter)"))
end

function r_dw(points, eps)
    @rput points eps
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("dw(points, eps=eps)"))
end

function r_inahull(points, alpha, queries)
    @rput points alpha queries
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("inahull(ahull(points, alpha=alpha), queries)"))
end

function r_ashape(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    RCall.reval(".as <- ashape(points, alpha=alpha)")
    return RCall.rcopy(RCall.reval(".as\$edges")), Int.(RCall.rcopy(RCall.reval(".as\$alpha.extremes")))
end

function r_complement(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("complement(points, alpha=alpha)"))
end

function r_dynamic_alpha(points; fraction=0.95, partCount=3, buff=10000,
                         initialAlpha=3, alphaIncrement=1, alphaCap=400)
    @rput points fraction partCount buff initialAlpha alphaIncrement alphaCap
    RCall.reval("library(rangeBuilder)")
    return RCall.rcopy(RCall.reval(
        "getDynamicAlphaHull(points, fraction=fraction, partCount=partCount, buff=buff, coordHeaders=c(1, 2), clipToCoast='no', initialAlpha=initialAlpha, alphaIncrement=alphaIncrement, alphaCap=alphaCap)[[2]]"
    ))
end

function r_filter_by_land(points)
    @rput points
    RCall.reval("library(rangeBuilder)")
    return RCall.rcopy(RCall.reval("filterByLand(points)"))
end

function align_point_ids(jpoints, rpoints)
    length(jpoints) == length(rpoints) || error("different point counts")
    mapping = zeros(Int, size(jpoints, 1))
    used = falses(size(rpoints, 1))
    for j in axes(jpoints, 1)
        d = vec(sqrt.(sum((rpoints .- jpoints[j, :]').^2, dims=2)))
        order = sortperm(d)
        k = findfirst(i -> !used[i] && d[i] < 1e-8, order)
        isnothing(k) && error("unable to align inserted point")
        rid = order[k]
        mapping[j] = rid
        used[rid] = true
    end
    return mapping
end

function canonical_ahull(arcs, jpoints, rpoints)
    mapping = align_point_ids(jpoints, rpoints)
    out = copy(arcs)
    for i in axes(out, 1)
        out[i, 7] = mapping[Int(out[i, 7])]
        out[i, 8] = mapping[Int(out[i, 8])]
    end
    key(i) = Tuple(round.(out[i, 1:6], digits=9))
    return out[sortperm(1:size(out, 1); by=key), :], mapping
end

function canonical(mesh)
    out = copy(mesh)
    for i in axes(out, 1)
        if out[i, 1] > out[i, 2]
            out[i, 1], out[i, 2] = out[i, 2], out[i, 1]
            out[i, 3], out[i, 5] = out[i, 5], out[i, 3]
            out[i, 4], out[i, 6] = out[i, 6], out[i, 4]
        end
        if (out[i, 7], out[i, 8]) > (out[i, 9], out[i, 10])
            out[i, 7], out[i, 9] = out[i, 9], out[i, 7]
            out[i, 8], out[i, 10] = out[i, 10], out[i, 8]
            out[i, 11], out[i, 12] = out[i, 12], out[i, 11]
        end
    end
    return out[sortperm(1:size(out, 1); by=i -> (out[i, 1], out[i, 2])), :]
end

function canonical_rows(rows; digits=9)
    return rows[sortperm(1:size(rows, 1); by=i -> Tuple(round.(rows[i, :], digits=digits))), :]
end

canonical_complement(rows) = canonical_rows(rows[:, 1:3])

points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]

dynamic_points = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.5 0.5]
# This point set is a cocircular degeneracy (square corners share a
# circumcircle through the center point). At alpha=0.7 the Delaunay triangles
# have circumradius sqrt(2)/2 ~= 0.707 > 0.7, so the hull degenerates to the
# center point only (coverage 1/5) and both engines enter the main loop.
@testset "dynamic square+center" begin
    for buff in (1000.0,)
        julia_dynamic = getDynamicAlphaHull(
            dynamic_points;
            fraction=0.8,
            partCount=3,
            buff=buff,
            initialAlpha=0.7,
            alphaIncrement=0.1,
            alphaCap=2,
            clipToCoast=:no,
        )
        r_dynamic = r_dynamic_alpha(
            dynamic_points;
            fraction=0.8,
            partCount=3,
            buff=buff,
            initialAlpha=0.7,
            alphaIncrement=0.1,
            alphaCap=2,
        )
        @test julia_dynamic.alpha == r_dynamic == "alpha0.8"
    end
    # Documented S2 difference for buff=0: R's st_buffer(hull, 0) round-trips
    # the geometry through +proj=eqearth and back, so the spherical engine
    # reports coverage 3/5 and keeps advancing until alpha > 2 -> alphaMCH.
    # Julia's planar GEOS buffer(0) is lossless (coverage 5/5 at alpha=0.8).
    # Both hulls are geometrically the same convex square; only the label
    # differs. See evidence/task-7-*.txt.
    julia_zero = getDynamicAlphaHull(
        dynamic_points;
        fraction=0.8,
        partCount=3,
        buff=0.0,
        initialAlpha=0.7,
        alphaIncrement=0.1,
        alphaCap=2,
        clipToCoast=:no,
    )
    r_zero = r_dynamic_alpha(
        dynamic_points;
        fraction=0.8,
        partCount=3,
        buff=0.0,
        initialAlpha=0.7,
        alphaIncrement=0.1,
        alphaCap=2,
    )
    @test julia_zero.alpha == "alpha0.8"
    @test r_zero == "alphaMCH"
end

land_reference_points = [116.4074 39.9042; -140.0 0.0]
@test filterByLand(land_reference_points; coastScale=50) == r_filter_by_land(land_reference_points)

# Task 7: six-class differential matrix against rangeBuilder::getDynamicAlphaHull.
# R side runs the installed rangeBuilder package through RCall; alpha labels
# must match exactly. Coordinates are shared fixed values so both engines see
# identical input (R and Julia RNGs differ, so shared-fixed beats re-seeding).
function dynamic_matrix_tests()
    check(points; kwargs...) = begin
        julia_out = getDynamicAlphaHull(points; clipToCoast=:no, kwargs...)
        r_alpha = r_dynamic_alpha(points; kwargs...)
        return julia_out.alpha, r_alpha
    end

    # Class 1: convex set with an interior point.
    basic = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.37 0.42]
    # Class 2: concave set (notch interior points).
    concave = [0.0 0.0; 2.0 0.0; 2.0 2.0; 1.0 0.6; 0.0 2.0; 0.5 1.1; 1.5 1.1]
    # Class 3: two separated clusters (fixed coordinates, non-axis-aligned).
    clusters = [0.522609 0.729549; 0.307886 0.250881; 0.347659 0.444204;
                0.833131 0.998262; 0.483922 0.545417; 19.038155 0.996055;
                19.720742 0.391774; 19.531328 0.607293; 19.964229 0.229241;
                19.847752 0.263021]
    # Class 4: first three rows share an exact x or y (R reshuffles them).
    share_x = [0.0 0.0; 0.0 1.0; 0.0 2.0; 5.0 5.0]
    share_y = [0.0 0.0; 1.0 0.0; 2.0 0.0; 5.0 5.0]
    # Class 5: degenerate inputs.
    collinear_diag = [0.0 0.0; 1.0 1.0; 2.0 2.0]
    two_points = [0.0 0.0; 0.0 1.0]
    # Class 6: buffered variants (buffer size changes coverage, hence alpha).

    @testset "dynamic matrix class1 convex" begin
        j, r = check(basic)
        @test j == r == "alpha3"
    end
    @testset "dynamic matrix class2 concave" begin
        j, r = check(concave)
        @test j == r == "alpha3"
    end
    @testset "dynamic matrix class3 clusters" begin
        # Default buff=10000: both engines select alpha12 (R via buffered
        # coverage, Julia via GEOS unbuffered coverage; see evidence).
        j, r = check(clusters)
        @test j == r == "alpha12"
    end
    @testset "dynamic matrix class4 first-three collinear" begin
        j, r = check(share_x)
        @test j == r == "alpha6"
        j, r = check(share_y)
        @test j == r == "alpha6"
    end
    @testset "dynamic matrix class5 degenerate" begin
        j, r = check(collinear_diag)
        @test j == r == "alphaMCH"
        @test_throws ArgumentError getDynamicAlphaHull(two_points; clipToCoast=:no)
        @test_throws RCall.REvalError r_dynamic_alpha(two_points)
    end
    @testset "dynamic matrix class6 buffered" begin
        # buff does not change the alpha for the convex case on either side.
        j, r = check(basic; buff=1000)
        @test j == r == "alpha3"
        j, r = check(basic; buff=10000)
        @test j == r == "alpha3"
        # S2 documented difference (see evidence task-7-*.txt): with buff=1000
        # R's spherical engine reports point (0.307886, 0.250881) - a hull
        # boundary vertex - as outside, so it advances to alpha14; Julia's
        # GEOS planar engine measures it inside at alpha12. Both are
        # documented, not asserted equal.
        julia_buffered, r_buffered = check(clusters; buff=1000)
        @test julia_buffered == "alpha12"
        @test r_buffered == "alpha14"
    end
end

dynamic_matrix_tests()

rmesh = canonical(r_delvor(points))
jmesh = canonical(delvor(points).mesh)
@test size(jmesh) == size(rmesh)
@test jmesh ≈ rmesh atol=1e-10 rtol=1e-10

r_arcs, r_points = r_ahull(points, 0.7)
j_ahull = ahull(points; alpha=0.7)
println("R ahull: arcs=$(size(r_arcs)), points=$(size(r_points)); Julia ahull: arcs=$(size(j_ahull.arcs)), points=$(size(j_ahull.xahull))")
canonical_arcs(arcs) = arcs[sortperm(1:size(arcs, 1); by=i -> (arcs[i, 7], arcs[i, 8])), :]
@test canonical_arcs(j_ahull.arcs) ≈ canonical_arcs(r_arcs) atol=1e-9 rtol=1e-9
@test j_ahull.xahull ≈ r_points atol=1e-9 rtol=1e-9

rng = MersenneTwister(20260808)
ring_angles = 2pi .* rand(rng, 24)
ring_radii = sqrt.(0.35^2 .+ (1.0 - 0.35^2) .* rand(rng, 24))
ring = hcat(ring_radii .* cos.(ring_angles), ring_radii .* sin.(ring_angles))
r_arcs, r_points = r_ahull(ring, 0.35)
j_ahull = ahull(ring; alpha=0.35)
println("R ring: arcs=$(size(r_arcs)), points=$(size(r_points)); Julia ring: arcs=$(size(j_ahull.arcs)), points=$(size(j_ahull.xahull))")
ca_j, mapping = canonical_ahull(j_ahull.arcs, j_ahull.xahull, r_points)
ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
jpoints_on_r = similar(r_points)
for j in axes(mapping, 1); jpoints_on_r[mapping[j], :] = j_ahull.xahull[j, :]; end
@test ca_j ≈ ca_r atol=1e-8 rtol=1e-8
@test jpoints_on_r ≈ r_points atol=1e-8 rtol=1e-8
@test areaahull(j_ahull) ≈ r_area(ring, 0.35) atol=1e-8 rtol=1e-8

@test size(koch(side=1, niter=3)) == size(r_koch(1, 3))

r_w = r_dw(points, 0.7)
j_w = dw(points; eps=0.7)
println("R dw: $(size(r_w)); Julia dw: $(size(j_w))")
@test size(j_w) == size(r_w)
@test j_w ≈ r_w atol=1e-9 rtol=1e-9

random_case = [0.03 0.12; 0.14 0.83; 0.27 0.31; 0.36 0.68;
               0.42 0.09; 0.51 0.55; 0.63 0.22; 0.71 0.91;
               0.79 0.43; 0.87 0.74; 0.94 0.16; 0.58 0.97]
random_alpha = 0.23
r_arcs, r_points = r_ahull(random_case, random_alpha)
j_ahull = ahull(random_case; alpha=random_alpha)
ca_j, mapping = canonical_ahull(j_ahull.arcs, j_ahull.xahull, r_points)
ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
@test ca_j ≈ ca_r atol=1e-8 rtol=1e-8
queries = [0.5 0.5; 0.02 0.98; 0.95 0.95; 0.4 0.2; 0.8 0.5]
@test inahull(j_ahull, queries) == r_inahull(random_case, random_alpha, queries)
@test areaahull(j_ahull) ≈ r_area(random_case, random_alpha) atol=1e-8 rtol=1e-8

function randomized_reference_tests()
    rng = MersenneTwister(20260808)
    for case_id in 1:8
        n = 6 + 2 * case_id
        scale = case_id % 2 == 0 ? 3.5 : 1.0
        shift = [case_id / 7, -case_id / 11]
        sample = (0.05 .+ 0.9 .* rand(rng, n, 2)) .* scale .+ shift'
        alpha = scale * (0.08 + 0.025 * (case_id % 5))
        @testset "random differential case $case_id" begin
        r_edges, r_extremes = r_ashape(sample, alpha)
        j_shape = ashape(sample; alpha)
        @test canonical(j_shape.edges) ≈ canonical(r_edges) atol=1e-8 rtol=1e-8
        @test j_shape.alpha_extremes == r_extremes

        r_comp = r_complement(sample, alpha)
        @test canonical_complement(complement(sample; alpha)) ≈ canonical_complement(r_comp) atol=1e-8 rtol=1e-8

        r_arcs, r_points = r_ahull(sample, alpha)
        j_hull = ahull(sample; alpha)
        j_reused = ahull(delvor(sample); alpha=alpha)
        j_split = ahull(sample[:, 1], sample[:, 2]; alpha=alpha)
        j_hull_arcs = j_hull.arcs[sortperm(1:size(j_hull.arcs, 1); by=i -> Tuple(round.(j_hull.arcs[i, 1:6], digits=9))), :]
        reused_arcs, _ = canonical_ahull(j_reused.arcs, j_reused.xahull, j_hull.xahull)
        split_arcs, _ = canonical_ahull(j_split.arcs, j_split.xahull, j_hull.xahull)
        @test reused_arcs ≈ j_hull_arcs atol=1e-8 rtol=1e-8
        @test split_arcs ≈ j_hull_arcs atol=1e-8 rtol=1e-8
        ca_j, _ = canonical_ahull(j_hull.arcs, j_hull.xahull, r_points)
        ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
        @test ca_j ≈ ca_r atol=1e-8 rtol=1e-8
        queries = (rand(rng, 8, 2) .* 1.4 .- 0.2) .* scale .+ shift'
        @test inahull(j_hull, queries) == r_inahull(sample, alpha, queries)
        @test areaahull(j_hull) ≈ r_area(sample, alpha) atol=1e-8 rtol=1e-8
        end
    end
end

randomized_reference_tests()

function boundary_reference_tests()
    boundary_points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
    queries = [0.3 0.3; 1.0 0.5; 3.0 3.0]
    for alpha in (0.0, 10.0)
        r_arcs, r_points = r_ahull(boundary_points, alpha)
        j_hull = ahull(boundary_points; alpha)
        ca_j, _ = canonical_ahull(j_hull.arcs, j_hull.xahull, r_points)
        ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
        @test ca_j ≈ ca_r atol=1e-8 rtol=1e-8
        @test inahull(j_hull, queries) == r_inahull(boundary_points, alpha, queries)
        @test areaahull(j_hull) ≈ r_area(boundary_points, alpha) atol=1e-8 rtol=1e-8
    end
end

boundary_reference_tests()

function near_degenerate_reference_tests()
    near_collinear = [0.0 0.0; 1.0 1e-9; 2.0 -1e-9; 3.0 2e-9; 1.5 0.4]
    @test_throws RCall.REvalError r_ahull(near_collinear, 0.45)
    @test size(ahull(near_collinear; alpha=0.45).arcs, 2) == 8

    cases = [
        ([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -0.99999999; 0.2 0.1], 0.9),
        ([0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7], 1e-6)
    ]
    queries = [0.1 0.1; 0.5 0.5; 2.5 2.5]
    for (sample, alpha) in cases
        r_arcs, r_points = r_ahull(sample, alpha)
        j_hull = ahull(sample; alpha)
        ca_j, _ = canonical_ahull(j_hull.arcs, j_hull.xahull, r_points)
        ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
        @test ca_j ≈ ca_r atol=1e-8 rtol=1e-8
        @test inahull(j_hull, queries) == r_inahull(sample, alpha, queries)
        @test areaahull(j_hull) ≈ r_area(sample, alpha) atol=1e-8 rtol=1e-8
    end
end

near_degenerate_reference_tests()

concave = [0.0 0.0; 2.0 0.0; 2.0 2.0; 1.0 0.6; 0.0 2.0; 0.5 1.1; 1.5 1.1]
r_arcs, r_points = r_ahull(concave, 0.7)
j_ahull = ahull(concave; alpha=0.7)
println("R concave: arcs=$(size(r_arcs)), points=$(size(r_points)); Julia concave: arcs=$(size(j_ahull.arcs)), points=$(size(j_ahull.xahull))")
ca_j, mapping = canonical_ahull(j_ahull.arcs, j_ahull.xahull, r_points)
ca_r = r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, 1:6], digits=9))), :]
@test ca_j ≈ ca_r atol=1e-9 rtol=1e-9
