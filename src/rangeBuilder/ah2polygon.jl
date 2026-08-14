"""Convert an `AHull` arc table to GeoInterface polygon geometry.

The alpha-hull core stores circular arcs in traversal order.  This adapter
samples each arc using the same density as R `ah2sf`, joins adjacent arcs by
their endpoint ids, and returns a polygon or multipolygon wrapper.
"""
function _ahull_arc_samples(row::AbstractVector; increment::Real, rnd::Integer)
    theta = row[6]
    n = max(2, 2 + round(Int, increment * (theta / 2)))
    pts = arc(row[1:2], row[3], row[4:5], theta; n)
    return round.(pts; digits=rnd)
end

function _orient_arc_samples(samples::AbstractMatrix, row::AbstractVector,
                             xahull::AbstractMatrix)
    a = Int(round(row[7]))
    b = Int(round(row[8]))
    if 1 <= a <= size(xahull, 1) && 1 <= b <= size(xahull, 1)
        forward = sum(abs2, samples[1, :] .- xahull[a, :]) +
                  sum(abs2, samples[end, :] .- xahull[b, :])
        reverse_cost = sum(abs2, samples[end, :] .- xahull[a, :]) +
                       sum(abs2, samples[1, :] .- xahull[b, :])
        reverse_cost < forward && return Base.reverse(samples; dims=1)
    end
    return samples
end

function _arc_components(h::AHull, rows::Vector{Int})
    # Equivalent to the endpoint-reordering state machine in R's ah2sf:
    # repeatedly attach an unused arc by its start endpoint, or by reversing
    # an arc whose end endpoint matches the current end.
    used = falses(size(h.arcs, 1))
    components = Vector{Vector{Tuple{Int,Bool}}}()
    for first_row in rows
        used[first_row] && continue
        used[first_row] = true
        component = Tuple{Int,Bool}[(first_row, false)]
        endpoint = Int(round(h.arcs[first_row, 8]))
        while true
            next_row = nothing
            reverse = false
            for candidate in rows
                used[candidate] && continue
                if Int(round(h.arcs[candidate, 7])) == endpoint
                    next_row = candidate
                    break
                end
            end
            if isnothing(next_row)
                for candidate in rows
                    used[candidate] && continue
                    if Int(round(h.arcs[candidate, 8])) == endpoint
                        next_row = candidate
                        reverse = true
                        break
                    end
                end
            end
            isnothing(next_row) && break
            row = Int(next_row)
            used[row] = true
            push!(component, (row, reverse))
            endpoint = reverse ? Int(round(h.arcs[row, 7])) : Int(round(h.arcs[row, 8]))
        end
        push!(components, component)
    end
    return components
end

function _join_arc_components(rows::Vector{Int}, h::AHull;
                              increment::Real, rnd::Integer, tol::Real)
    components = Vector{Vector{Tuple{Float64,Float64}}}()
    isempty(rows) && return components
    # R ah2sf reorders the arc table with a state machine, then joins arcs in
    # one pass over that reordered sequence.  `_arc_components` reproduces the
    # reordering; flattening its component orders yields the same sequence.
    order = Tuple{Int,Bool}[]
    for component_order in _arc_components(h, rows)
        append!(order, component_order)
    end
    # Join arcs by coordinate continuity, exactly as R does:
    #  * an arc whose first sample is continuous with the current line's last
    #    sample is appended (minus its duplicated first point);
    #  * a discontinuous arc closes the current line and is itself dropped
    #    (its samples never enter a line), matching R's `prevx <- NULL`;
    #  * a line is closed by replacing its last sample with its first sample,
    #    as R does with `prevx[length(prevx)] <- prevx[1]`; the final arc in
    #    the sequence closes only when it is continuous.
    current = Tuple{Float64,Float64}[]
    n = length(order)
    for (position, (i, flipped)) in enumerate(order)
        row = copy(h.arcs[i, :])
        if flipped
            row[7], row[8] = row[8], row[7]
        end
        samples = _orient_arc_samples(_ahull_arc_samples(row; increment, rnd), row, h.xahull)
        points = [(Float64(p[1]), Float64(p[2])) for p in eachrow(samples)]
        if isempty(current)
            append!(current, points)
        elseif _arc_samples_continuous(points[1], current[end]; rnd, tol)
            append!(current, points[2:end])
            if position == n
                current[end] = current[1]
                push!(components, current)
            end
        else
            current[end] = current[1]
            push!(components, current)
            current = Tuple{Float64,Float64}[]
        end
    end
    # R's badLines step: after joining, any line with fewer than four points
    # (sf::st_cast(line, "POINT") followed by nrow < 4) is removed.  This is
    # what drops three-point degenerate rings that a pure endpoint-id join
    # would otherwise retain.
    filter!(line -> length(line) >= 4, components)
    return components
end

function _arc_samples_continuous(first::Tuple{Float64,Float64},
                                 last::Tuple{Float64,Float64};
                                 rnd::Integer, tol::Real)
    (fx, fy) = first
    (lx, ly) = last
    x_ok = fx == round(lx; digits=rnd) || abs(fx - lx) < tol
    y_ok = fy == round(ly; digits=rnd) || abs(fy - ly) < tol
    return x_ok && y_ok
end

"""Convert an alpha hull to a GeoInterface `Polygon` or `MultiPolygon`.

`increment` and `rnd` mirror the corresponding R `ah2sf` arguments.  Empty
or point-only alpha hulls return `nothing`, matching R's empty geometry result.
"""
function ah2polygon(h::AHull; increment::Real=360, rnd::Integer=10,
                    crs=4326, tol::Real=1e-4)
    increment > 0 || throw(ArgumentError("increment must be positive"))
    rnd >= 0 || throw(ArgumentError("rnd must be nonnegative"))
    tol >= 0 || throw(ArgumentError("tol must be nonnegative"))
    rows = findall(>(0), h.arcs[:, 3])
    components = _join_arc_components(rows, h; increment, rnd, tol)
    isempty(components) && return nothing
    wrappers = GeoInterface.Wrappers
    polygons = [wrappers.Polygon([ring]; crs) for ring in components]
    return length(polygons) == 1 ? polygons[1] : wrappers.MultiPolygon(polygons; crs)
end
