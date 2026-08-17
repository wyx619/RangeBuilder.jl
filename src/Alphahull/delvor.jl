struct DelVor
    x::Matrix{Float64}
    mesh::Matrix{Float64}
    triangulation
end

function _points(x, y=nothing)
    if y === nothing
        a = Matrix{Float64}(x)
        size(a, 2) == 2 || throw(ArgumentError("x must have two columns"))
    else
        length(x) == length(y) || throw(ArgumentError("x and y must have equal length"))
        a = hcat(Float64.(x), Float64.(y))
    end
    size(a, 1) >= 3 || throw(ArgumentError("at least three points are required"))
    any(v -> !isfinite(v), a) && throw(ArgumentError("coordinates must be finite"))
    seen = Set{NTuple{2,Float64}}()
    sizehint!(seen, size(a, 1))
    for i in axes(a, 1)
        point = (a[i, 1], a[i, 2])
        point in seen && throw(ArgumentError("duplicate data points"))
        push!(seen, point)
    end
    base = a[1, :]
    ref = 0
    for i in 2:size(a, 1)
        if a[i, 1] != base[1] || a[i, 2] != base[2]
            ref = i
            break
        end
    end
    ref == 0 && throw(ArgumentError("points must not be collinear"))
    cross = (a[ref, 1] - base[1]) .* (a[:, 2] .- base[2]) .-
            (a[ref, 2] - base[2]) .* (a[:, 1] .- base[1])
    maximum(abs, cross) > eps(Float64) || throw(ArgumentError("points must not be collinear"))
    return a
end

function _dummycoor(x, i, j, m, away, interior)
    vx, vy = x[j,2] - x[i,2], -(x[j,1] - x[i,1])
    norm2 = vx^2 + vy^2
    norm2 > 0 && ((vx, vy) = (vx / norm2, vy / norm2))
    side = (x[j,1] - x[i,1]) * (x[interior,2] - x[i,2]) -
           (x[j,2] - x[i,2]) * (x[interior,1] - x[i,1])
    return side < 0 ? (m[1] - away * vx, m[2] - away * vy) :
                      (m[1] + away * vx, m[2] + away * vy)
end

"""Construct the Delaunay/Voronoi edge table used by `ashape`.

`backend=:delaunay` retains the existing high-performance JuliaGeometry
implementation. `backend=:shull` is an experimental, pure-Julia port of the
SHull path used by R `interp::tri.mesh()`; it is provided for differential
validation and does not require R or RCall at runtime.
"""
function delvor(x, y=nothing; rng=Random.MersenneTwister(0),
                randomise::Bool=false, backend::Symbol=:delaunay)
    xy = _points(x, y)
    if backend === :shull
        randomise && throw(ArgumentError("the deterministic SHull backend does not support randomise=true"))
        isnothing(_shull_float32_duplicate_pair(xy)) ||
            throw(ArgumentError("duplicate data points"))
        shull = _shull_final_triads(xy; rng=rng)
        return DelVor(xy, _shull_mesh(xy, shull), shull)
    end
    backend === :delaunay || throw(ArgumentError("backend must be :delaunay or :shull"))
    pts = [(xy[i, 1], xy[i, 2]) for i in axes(xy, 1)]
    # R's interp::tri.mesh builds a deterministic triangulation. Randomised
    # insertion can rotate otherwise identical alpha-hull arc cycles; because
    # rangeBuilder::ah2sf is row-order sensitive, that changes its polygons.
    # `randomise=true` retains the previous Julia traversal option for callers
    # who explicitly need it; R-compatible construction is the default.
    tri = triangulate(pts; rng, randomise)
    span = max(maximum(xy[:, 1]) - minimum(xy[:, 1]), maximum(xy[:, 2]) - minimum(xy[:, 2]))
    span = span > 0 ? span : 1.0
    edges = sort!(collect(DelaunayTriangulation.each_solid_edge(tri)))
    rows = Matrix{Float64}(undef, length(edges), 12)
    for (row, e) in enumerate(edges)
        i, j = e
        k1 = DelaunayTriangulation.get_adjacent(tri, i, j)
        k2 = DelaunayTriangulation.get_adjacent(tri, j, i)
        if k1 > 0
            c1 = DelaunayTriangulation.triangle_circumcenter(tri, (i, j, k1))
        else
            c1 = DelaunayTriangulation.triangle_circumcenter(tri, (j, i, k2))
        end
        if k1 > 0 && k2 > 0
            c2 = DelaunayTriangulation.triangle_circumcenter(tri, (j, i, k2)); bp1 = 0.0; bp2 = 0.0
        else
            interior = k1 > 0 ? k1 : k2
            c2 = _dummycoor(xy, i, j, c1, span, interior); bp1 = 0.0; bp2 = 1.0
        end
        @inbounds begin
            rows[row, 1] = i
            rows[row, 2] = j
            rows[row, 3] = pts[i][1]
            rows[row, 4] = pts[i][2]
            rows[row, 5] = pts[j][1]
            rows[row, 6] = pts[j][2]
            rows[row, 7] = c1[1]
            rows[row, 8] = c1[2]
            rows[row, 9] = c2[1]
            rows[row, 10] = c2[2]
            rows[row, 11] = bp1
            rows[row, 12] = bp2
        end
    end
    return DelVor(xy, rows, tri)
end
