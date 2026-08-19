# Optional development-only differential tests against the locally installed,
# unmodified R packages `alphahull` and `rangeBuilder`.
#
# Run explicitly with: julia --project=test test/reference_r.jl
# This file is deliberately not included by `test/runtests.jl`: RCall and R are
# development dependencies, not runtime requirements of RangeBuilder.jl.
using Test
using Random
using RCall
using RangeBuilder

RCall.reval("options(warn = -1)")

function r_delvor(points)
    @rput points
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("delvor(points)\$mesh"))
end

function r_ahull(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    RCall.reval(".reference_hull <- ahull(points, alpha=alpha)")
    arcs = RCall.rcopy(RCall.reval(".reference_hull\$arcs"))
    xahull = RCall.rcopy(RCall.reval(".reference_hull\$xahull"))
    return arcs, xahull
end

function r_area(points, alpha)
    @rput points alpha
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("areaahull(ahull(points, alpha=alpha))"))
end

function r_inahull(points, alpha, queries)
    @rput points alpha queries
    RCall.reval("library(alphahull)")
    return RCall.rcopy(RCall.reval("inahull(ahull(points, alpha=alpha), queries)"))
end

function r_dynamic_alpha(points; fraction=0.95, partCount=3, buff=10000,
                         initialAlpha=3, alphaIncrement=1, alphaCap=400)
    @rput points fraction partCount buff initialAlpha alphaIncrement alphaCap
    RCall.reval("library(rangeBuilder)")
    return RCall.rcopy(RCall.reval(
        "getDynamicAlphaHull(points, fraction=fraction, partCount=partCount, " *
        "buff=buff, coordHeaders=c(1, 2), clipToCoast='no', " *
        "initialAlpha=initialAlpha, alphaIncrement=alphaIncrement, " *
        "alphaCap=alphaCap)[[2]]"
    ))
end

function canonical_mesh(mesh)
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
    return out[sortperm(1:size(out, 1); by=i -> Tuple(round.(out[i, :], digits=9))), :]
end

function align_point_ids(jpoints, rpoints)
    length(jpoints) == length(rpoints) || error("different point counts")
    mapping = zeros(Int, size(jpoints, 1))
    used = falses(size(rpoints, 1))
    for j in axes(jpoints, 1)
        distances = vec(sqrt.(sum((rpoints .- jpoints[j, :]').^2, dims=2)))
        order = sortperm(distances)
        k = findfirst(i -> !used[i] && distances[i] < 1e-8, order)
        isnothing(k) && error("unable to align inserted point")
        mapping[j] = order[k]
        used[order[k]] = true
    end
    return mapping
end

function canonical_arcs(arcs, jpoints, rpoints)
    mapping = align_point_ids(jpoints, rpoints)
    out = copy(arcs)
    for i in axes(out, 1)
        out[i, 7] = mapping[Int(out[i, 7])]
        out[i, 8] = mapping[Int(out[i, 8])]
    end
    return out[sortperm(1:size(out, 1); by=i -> Tuple(round.(out[i, :], digits=9))), :]
end

function assert_reference_equivalence(points, alpha, queries)
    rmesh = canonical_mesh(r_delvor(points))
    jmesh = canonical_mesh(delvor(points).mesh)
    # The mesh is an implementation detail. Different valid Delaunay
    # diagonals can yield different rows even before alpha-hull construction.
    @test size(jmesh) == size(rmesh)

    r_arcs, r_points = r_ahull(points, alpha)
    j_hull = ahull(points; alpha)
    @test canonical_arcs(j_hull.arcs, j_hull.xahull, r_points) ≈
          r_arcs[sortperm(1:size(r_arcs, 1); by=i -> Tuple(round.(r_arcs[i, :], digits=9))), :] atol=1e-8 rtol=1e-8
    @test inahull(j_hull, queries) == r_inahull(points, alpha, queries)
    @test areaahull(j_hull) ≈ r_area(points, alpha) atol=1e-8 rtol=1e-8
end

@testset "Original R reference: nondegenerate geometry" begin
    basic = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
    assert_reference_equivalence(basic, 0.7, [0.3 0.3; 1.0 0.5; 3.0 3.0])

    concave = [0.0 0.0; 2.0 0.0; 2.0 2.0; 1.0 0.6; 0.0 2.0; 0.5 1.1; 1.5 1.1]
    assert_reference_equivalence(concave, 0.7, [0.3 0.3; 1.0 1.4; 3.0 3.0])

    rng = MersenneTwister(20260808)
    points = 0.05 .+ 0.9 .* rand(rng, 20, 2)
    assert_reference_equivalence(points, 0.17, rand(rng, 8, 2) .* 1.4 .- 0.2)
end

@testset "Original R reference: stable dynamic alpha" begin
    points = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.37 0.42]
    julia_result = getDynamicAlphaHull(points; clipToCoast=:no)
    r_result = r_dynamic_alpha(points)
    @test julia_result.alpha == r_result == "alpha3"
end
