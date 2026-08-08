struct CircleIntersection
    n_cut::Float64
    v1::NTuple{2,Float64}
    theta1::Float64
    v2::NTuple{2,Float64}
    theta2::Float64
end

"""Port of R `alphahull::inter` for two circles."""
function inter(c11::Real, c12::Real, r1::Real, c21::Real, c22::Real, r2::Real)
    d = hypot(c21 - c11, c22 - c12)
    if d == 0 && r1 == r2
        return CircleIntersection(Inf, (0.0, 0.0), 0.0, (0.0, 0.0), 0.0)
    elseif d > r1 + r2 || d < max(r1, r2) - min(r1, r2)
        return CircleIntersection(0.0, (0.0, 0.0), 0.0, (0.0, 0.0), 0.0)
    elseif d == max(r1, r2) - min(r1, r2) || d == r1 + r2
        return CircleIntersection(1.0, (0.0, 0.0), 0.0, (0.0, 0.0), 0.0)
    end
    d1 = (d^2 - r2^2 + r1^2) / (2d)
    d2 = d - d1
    vx, vy = (c21 - c11) / d, (c22 - c12) / d
    return CircleIntersection(2.0, (vx, vy), acos(d1 / r1), (-vx, -vy), acos(d2 / r2))
end
