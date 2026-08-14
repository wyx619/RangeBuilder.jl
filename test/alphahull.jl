@testset "Alpha-hull core" begin
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
    @test areaahull(ah) >= 0

    ci = inter(0, 0, 1, 1, 0, 1)
    @test ci.n_cut == 2
    @test ci.v1 == (1.0, 0.0)
    @test ci.theta1 ≈ pi / 3
    @test inter(0, 0, 1, 3, 0, 1).n_cut == 0
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
end
