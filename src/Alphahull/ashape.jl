struct AShape
    edges::Matrix{Float64}
    length::Float64
    alpha::Float64
    alpha_extremes::Vector{Int}
    delvor::DelVor
    x::Matrix{Float64}
end

# R's `grDevices::chull` orders convex-hull vertices by the angle about the
# hull centroid. DelaunayTriangulation supplies the hull itself; this only
# reproduces R's public ordering because `dw` exposes the alpha extremes.
function _r_hull_indices(x, tri)
    hull = reverse(unique(Int.(DelaunayTriangulation.get_convex_hull_vertices(tri))))
    center_x = sum(@view x[hull, 1]) / length(hull)
    center_y = sum(@view x[hull, 2]) / length(hull)
    sort!(hull; by=i -> atan(x[i, 2] - center_y, center_x - x[i, 1]))
    return hull
end

function _r_hull_indices(x, shull::NamedTuple{(:points, :triads, :hull)})
    hull = [point.id for point in shull.hull]
    center_x = sum(@view x[hull, 1]) / length(hull)
    center_y = sum(@view x[hull, 2]) / length(hull)
    sort!(hull; by=i -> atan(x[i, 2] - center_y, center_x - x[i, 1]))
    return hull
end

"""Compute the alpha-shape using the same edge criterion as R `alphahull`."""
function _ashape(x, y, alpha::Real)
    alpha >= 0 || throw(ArgumentError("alpha must be nonnegative"))
    dv = x isa DelVor && y === nothing ? x : delvor(x, y)
    m, xy = dv.mesh, dv.x
    n = size(xy, 1)
    dm1 = hypot.(m[:,3] .- m[:,7], m[:,4] .- m[:,8])
    dm2 = hypot.(m[:,3] .- m[:,9], m[:,4] .- m[:,10])
    dm1[m[:,11] .== 1] .= Inf; dm2[m[:,12] .== 1] .= Inf
    hull = _r_hull_indices(xy, dv.triangulation)
    ext = Set(hull)
    maxd = zeros(Float64, n)
    for r in axes(m,1)
        i, j = Int(m[r,1]), Int(m[r,2])
        maxd[i] = max(maxd[i], dm1[r], dm2[r])
        maxd[j] = max(maxd[j], dm1[r], dm2[r])
    end
    interior = [i for i in 1:n if !(i in ext)]
    append!(hull, [i for i in interior if maxd[i] > alpha])
    union!(ext, hull)
    keep = BitVector(undef, size(m,1))
    for r in axes(m,1)
        i, j = Int(m[r,1]), Int(m[r,2])
        pmx, pmy = (m[r,3] + m[r,5]) / 2, (m[r,4] + m[r,6]) / 2
        between = ((m[r,7] == m[r,9] && min(m[r,8],m[r,10]) <= pmy <= max(m[r,8],m[r,10])) ||
                   (m[r,7] != m[r,9] && min(m[r,7],m[r,9]) <= pmx <= max(m[r,7],m[r,9])))
        half = hypot(m[r,3] - m[r,5], m[r,4] - m[r,6]) / 2
        lmin = min(dm1[r], dm2[r], between ? half : Inf)
        lmax = max(dm1[r], dm2[r], between ? half : -Inf)
        keep[r] = (i in ext) && (j in ext) && lmin <= alpha <= lmax
    end
    edges = m[keep, :]
    len = sum(hypot.(edges[:,3] .- edges[:,5], edges[:,4] .- edges[:,6]))
    return AShape(edges, len, Float64(alpha), hull, dv, xy)
end

ashape(x, y, alpha::Real) = _ashape(x, y, alpha)
ashape(x; alpha::Real) = _ashape(x, nothing, alpha)
ashape(x, y; alpha::Real) = _ashape(x, y, alpha)
