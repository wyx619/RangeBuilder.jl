using DelimitedFiles

const GOLDEN_DIR = joinpath(@__DIR__, "fixtures", "generated")

function _golden_csv(tag, suffix)
    path = joinpath(GOLDEN_DIR, "$(tag)_$(suffix).csv")
    length(readlines(path)) == 1 && return zeros(Float64, 0, suffix == "ashape_edges" ? 12 : suffix == "complement" ? 19 : 8)
    return Float64.(readdlm(path, ',', skipstart=1))
end

function _canonical_mesh(mesh)
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

function _rowset_approx(a, b; atol, rtol)
    size(a) == size(b) || return false
    used = falses(size(b, 1))
    for i in axes(a, 1)
        found = false
        for j in axes(b, 1)
            if !used[j] && all(isapprox(a[i, k], b[j, k]; atol, rtol) for k in axes(a, 2))
                used[j] = true
                found = true
                break
            end
        end
        found || return false
    end
    return true
end

function _align_golden_points(jpoints, rpoints)
    mapping = zeros(Int, size(jpoints, 1))
    used = falses(size(rpoints, 1))
    for j in axes(jpoints, 1)
        d = vec(sqrt.(sum((rpoints .- jpoints[j, :]').^2, dims=2)))
        for rid in sortperm(d)
            if !used[rid] && d[rid] < 1e-8
                mapping[j] = rid
                used[rid] = true
                break
            end
        end
        mapping[j] > 0 || error("unable to align golden alpha-hull point")
    end
    return mapping
end

function _canonical_arcs(arcs, jpoints, rpoints)
    mapping = _align_golden_points(jpoints, rpoints)
    out = copy(arcs)
    for i in axes(out, 1)
        out[i, 7] = mapping[Int(out[i, 7])]
        out[i, 8] = mapping[Int(out[i, 8])]
    end
    return out[sortperm(1:size(out, 1); by=i -> Tuple(round.(out[i, 1:6], digits=9))), :]
end

function run_golden_fixtures()
for (name, alpha) in (("square_center", 0.35), ("square_center", 0.70),
                      ("concave", 0.35), ("concave", 0.70),
                      ("ring", 0.35), ("ring", 0.70))
    alpha_tag = alpha == 0.35 ? "0.35" : "0.70"
    tag = "$(name)_alpha_$(alpha_tag)"
    points = _golden_csv(tag, "points")
    rmesh = _golden_csv(tag, "mesh")
    redges = _golden_csv(tag, "ashape_edges")
    rcomp = _golden_csv(tag, "complement")
    rarcs = _golden_csv(tag, "ahull_arcs")
    rpoints = _golden_csv(tag, "ahull_points")
    rsummary = _golden_csv(tag, "summary")

    dv = delvor(points)
    @test size(dv.mesh) == size(rmesh)
    @test all(isfinite, dv.mesh)
    @test all(1 .<= dv.mesh[:, 1:2] .<= size(points, 1))

    shape = ashape(dv; alpha)
    @test _canonical_mesh(shape.edges) ≈ _canonical_mesh(redges) atol=1e-9 rtol=1e-9

    comp = complement(dv; alpha)
    @test size(comp) == size(rcomp)
    @test _rowset_approx(comp[:, 1:3], rcomp[:, 1:3]; atol=1e-8, rtol=1e-8)

    hull = ahull(dv; alpha)
    @test _canonical_arcs(hull.arcs, hull.xahull, rpoints) ≈
          _canonical_arcs(rarcs, rpoints, rpoints) atol=1e-8 rtol=1e-8
    @test _rowset_approx(hull.xahull, rpoints; atol=1e-9, rtol=1e-9)
    @test hull.length ≈ rsummary[1] atol=1e-9 rtol=1e-9
end
end

run_golden_fixtures()
