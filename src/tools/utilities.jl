"""Clockwise rotation used throughout the R package."""
rotation(v::AbstractVector, theta::Real) = [cos(theta) * v[1] + sin(theta) * v[2],
                                             -sin(theta) * v[1] + cos(theta) * v[2]]

function anglesArc(v::AbstractVector, theta::Real)
    thetaox = v[2] >= 0 ? acos(clamp(v[1], -1.0, 1.0)) : 2pi - acos(clamp(v[1], -1.0, 1.0))
    return [thetaox - theta, thetaox + theta]
end

function koch(; side=3.0, niter=5)
    niter >= 1 || throw(ArgumentError("niter must be at least one"))
    npoints = Int(3 * 4^(niter - 1))
    # The R implementation constructs the vertices in-place; this equivalent
    # refinement keeps the same clockwise equilateral seed and point order.
    p = [(-side / 2, 0.0), (side / 2, 0.0),
         (0.0, side * sqrt(3) / 2)]
    for _ in 2:niter
        q = NTuple{2,Float64}[]
        for i in eachindex(p)
            a, b = p[i], p[mod1(i + 1, 3)]
            v = ((b[1]-a[1])/3, (b[2]-a[2])/3)
            c = (a[1]+v[1], a[2]+v[2])
            d = (a[1]+2v[1], a[2]+2v[2])
            peak = (c[1] + cos(-pi/3)*v[1] + sin(-pi/3)*v[2],
                    c[2] - sin(-pi/3)*v[1] + cos(-pi/3)*v[2])
            append!(q, (a, c, peak, d))
        end
        p = q
    end
    return reduce(vcat, ([x y] for (x,y) in p))
end

function _inside_polygon(x, y, poly)
    inside = false
    j = size(poly, 1)
    for i in 1:size(poly, 1)
        xi, yi = poly[i,1], poly[i,2]
        xj, yj = poly[j,1], poly[j,2]
        ((yi > y) != (yj > y) && x < (xj-xi)*(y-yi)/(yj-yi)+xi) && (inside = !inside)
        j = i
    end
    return inside
end

function rkoch(n::Integer; side=3.0, niter=5, rng=Random.default_rng())
    n >= 0 || throw(ArgumentError("n must be nonnegative"))
    poly = koch(side=side, niter=niter)
    xmin, xmax = extrema(poly[:,1]); ymin, ymax = extrema(poly[:,2])
    out = Matrix{Float64}(undef, n, 2)
    count = 0
    while count < n
        x, y = rand(rng) * (xmax-xmin) + xmin, rand(rng) * (ymax-ymin) + ymin
        _inside_polygon(x, y, poly) || continue
        count += 1; out[count,:] .= (x, y)
    end
    return out
end
