"""Test whether each query point belongs to an alpha-convex hull."""
function inahull(h::AHull, points)
    p = points isa AbstractVector && length(points) == 2 ? reshape(Float64.(points), 1, 2) : Matrix{Float64}(points)
    size(p, 2) == 2 || throw(ArgumentError("points must have two columns"))
    outside = falses(size(p, 1))
    cp = h.complement
    for row in axes(cp, 1)
        r = cp[row, 3]
        if r < 0
            a, b, sig = cp[row, 1], cp[row, 2], Int(cp[row, 3])
            if sig == -4
                outside .|= p[:, 1] .< a
            elseif sig == -3
                outside .|= p[:, 1] .> a
            elseif sig == -2
                outside .|= p[:, 2] .< a .+ b .* p[:, 1]
            elseif sig == -1
                outside .|= p[:, 2] .> a .+ b .* p[:, 1]
            end
        elseif r > 0
            outside .|= hypot.(p[:, 1] .- cp[row, 1], p[:, 2] .- cp[row, 2]) .< r
        end
    end
    return .!outside
end
