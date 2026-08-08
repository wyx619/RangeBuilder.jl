using Test
using AlphaHull
using RCall
using DelimitedFiles
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

ring_path = joinpath(@__DIR__, "fixtures", "generated", "ring_alpha_0.35_points.csv")
ring = Float64.(readdlm(ring_path, ',', skipstart=1))
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
