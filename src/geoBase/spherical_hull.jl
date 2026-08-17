"""Return the counter-clockwise spherical convex-hull vertex indices.

The input is longitude/latitude in degrees. Points are converted with
`GeometryOps.UnitSpherical`, then gnomonically projected onto a tangent plane.
Great circles become lines under this projection, so a planar monotone-chain
hull has the same vertices as the spherical hull while all points lie in one
open hemisphere.
"""
function _spherical_convex_hull_indices(points::AbstractMatrix)
    count = size(points, 1)
    count >= 3 || throw(ArgumentError("at least three points are required"))

    to_sphere = GeometryOps.UnitSpherical.UnitSphereFromGeographic()
    unit_points = Vector{NTuple{3,Float64}}(undef, count)
    cx = 0.0
    cy = 0.0
    cz = 0.0
    @inbounds for index in 1:count
        point = to_sphere((Float64(points[index, 1]), Float64(points[index, 2])))
        value = (Float64(point.x), Float64(point.y), Float64(point.z))
        unit_points[index] = value
        cx += value[1]
        cy += value[2]
        cz += value[3]
    end

    norm_center = sqrt(cx * cx + cy * cy + cz * cz)
    norm_center > 0 || throw(ArgumentError("spherical hull has no unique hemisphere"))
    center = (cx / norm_center, cy / norm_center, cz / norm_center)
    all(point -> point[1] * center[1] + point[2] * center[2] + point[3] * center[3] > 64eps(Float64), unit_points) ||
        throw(ArgumentError("spherical hull points do not lie in one open hemisphere"))

    reference = abs(center[3]) < 0.9 ? (0.0, 0.0, 1.0) : (0.0, 1.0, 0.0)
    e1x = reference[2] * center[3] - reference[3] * center[2]
    e1y = reference[3] * center[1] - reference[1] * center[3]
    e1z = reference[1] * center[2] - reference[2] * center[1]
    norm_e1 = sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
    e1 = (e1x / norm_e1, e1y / norm_e1, e1z / norm_e1)
    e2 = (
        center[2] * e1[3] - center[3] * e1[2],
        center[3] * e1[1] - center[1] * e1[3],
        center[1] * e1[2] - center[2] * e1[1],
    )

    projected = Vector{NTuple{3,Float64}}(undef, count)
    @inbounds for index in 1:count
        point = unit_points[index]
        denominator = point[1] * center[1] + point[2] * center[2] + point[3] * center[3]
        projected[index] = (
            (point[1] * e1[1] + point[2] * e1[2] + point[3] * e1[3]) / denominator,
            (point[1] * e2[1] + point[2] * e2[2] + point[3] * e2[3]) / denominator,
            index,
        )
    end
    sort!(projected; by=point -> (point[1], point[2], point[3]))

    cross(origin, first, second) =
        (first[1] - origin[1]) * (second[2] - origin[2]) -
        (first[2] - origin[2]) * (second[1] - origin[1])
    lower = NTuple{3,Float64}[]
    for point in projected
        while length(lower) >= 2 && cross(lower[end - 1], lower[end], point) <= 0
            pop!(lower)
        end
        push!(lower, point)
    end
    upper = NTuple{3,Float64}[]
    for point in Iterators.reverse(projected)
        while length(upper) >= 2 && cross(upper[end - 1], upper[end], point) <= 0
            pop!(upper)
        end
        push!(upper, point)
    end
    hull = vcat(lower[1:end - 1], upper[1:end - 1])
    length(hull) >= 3 || throw(ArgumentError("spherical hull points are collinear"))
    return Int[point[3] for point in hull]
end

function _spherical_convex_hull_polygon(points::AbstractMatrix; crs=4326)
    indices = _spherical_convex_hull_indices(points)
    ring = [(Float64(points[index, 1]), Float64(points[index, 2])) for index in indices]
    push!(ring, first(ring))
    return GeoInterface.Wrappers.Polygon([ring]; crs)
end
