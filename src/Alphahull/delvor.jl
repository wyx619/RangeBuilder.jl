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

"""Construct R-compatible SHull Delaunay/Voronoi edge data used by `ashape`.

This is a pure-Julia port of the SHull path used by R `interp::tri.mesh()`.
It does not require R or RCall at runtime.
"""
function delvor(x, y=nothing; rng=Random.MersenneTwister(0),
                randomise::Bool=false)
    xy = _points(x, y)
    randomise && throw(ArgumentError("the deterministic SHull backend does not support randomise=true"))
    isnothing(_shull_float32_duplicate_pair(xy)) ||
        throw(ArgumentError("duplicate data points"))
    shull = _shull_final_triads(xy; rng=rng)
    return DelVor(xy, _shull_mesh(xy, shull), shull)
end
