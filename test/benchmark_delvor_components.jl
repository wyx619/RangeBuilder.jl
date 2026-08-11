using RangeBuilder
using DelaunayTriangulation
using Random

function measure(label, f; repeats=3)
    f()
    GC.gc()
    samples = [@timed f() for _ in 1:repeats]
    order = sortperm(getfield.(samples, :time))
    sample = samples[order[cld(repeats, 2)]]
    println(rpad(label, 28), " ", round(sample.time; digits=4), " s   ",
            round(sample.bytes / 2.0^20; digits=2), " MiB")
    return sample.value
end

function mesh_tail(xy, tri)
    pts = [(xy[i, 1], xy[i, 2]) for i in axes(xy, 1)]
    span = max(maximum(xy[:, 1]) - minimum(xy[:, 1]), maximum(xy[:, 2]) - minimum(xy[:, 2]))
    span = span > 0 ? span : 1.0
    edges = sort!(collect(DelaunayTriangulation.each_solid_edge(tri)))
    rows = Matrix{Float64}(undef, length(edges), 12)
    for (row, e) in enumerate(edges)
        i, j = e
        k1 = DelaunayTriangulation.get_adjacent(tri, i, j)
        k2 = DelaunayTriangulation.get_adjacent(tri, j, i)
        c1 = k1 > 0 ? DelaunayTriangulation.triangle_circumcenter(tri, (i, j, k1)) :
             DelaunayTriangulation.triangle_circumcenter(tri, (j, i, k2))
        if k1 > 0 && k2 > 0
            c2 = DelaunayTriangulation.triangle_circumcenter(tri, (j, i, k2))
        else
            interior = k1 > 0 ? k1 : k2
            vx, vy = xy[j, 2] - xy[i, 2], -(xy[j, 1] - xy[i, 1])
            norm2 = vx^2 + vy^2
            norm2 > 0 && ((vx, vy) = (vx / norm2, vy / norm2))
            side = (xy[j, 1] - xy[i, 1]) * (xy[interior, 2] - xy[i, 2]) -
                   (xy[j, 2] - xy[i, 2]) * (xy[interior, 1] - xy[i, 1])
            c2 = side < 0 ? (c1[1] - span * vx, c1[2] - span * vy) :
                             (c1[1] + span * vx, c1[2] + span * vy)
        end
        @inbounds begin
            rows[row, 1] = i; rows[row, 2] = j
            rows[row, 3] = pts[i][1]; rows[row, 4] = pts[i][2]
            rows[row, 5] = pts[j][1]; rows[row, 6] = pts[j][2]
            rows[row, 7] = c1[1]; rows[row, 8] = c1[2]
            rows[row, 9] = c2[1]; rows[row, 10] = c2[2]
            rows[row, 11] = 0.0; rows[row, 12] = (k1 > 0 && k2 > 0) ? 0.0 : 1.0
        end
    end
    return rows
end

for n in (4_000, 8_000, 10_000, 20_000)
    rng = MersenneTwister(20260809 + n)
    xy = rand(rng, n, 2)
    pts = [(xy[i, 1], xy[i, 2]) for i in axes(xy, 1)]
    println("\nn = ", n)
    tri = measure("triangulate (randomised)", () -> triangulate(pts; rng=MersenneTwister(7)))
    measure("mesh construction tail", () -> mesh_tail(xy, tri))
    measure("triangulate (input order)", () -> triangulate(pts; randomise=false))
end
