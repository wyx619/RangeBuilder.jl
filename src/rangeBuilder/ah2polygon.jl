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
    for component_order in _arc_components(h, rows)
        current = Tuple{Float64,Float64}[]
        previous_end = 0
        for (position_in_component, (i, flipped)) in enumerate(component_order)
            row = copy(h.arcs[i, :])
            if flipped
                row[7], row[8] = row[8], row[7]
            end
            samples = _orient_arc_samples(_ahull_arc_samples(row; increment, rnd), row, h.xahull)
            points = [(Float64(p[1]), Float64(p[2])) for p in eachrow(samples)]
            if position_in_component == 1 || Int(round(row[7])) == previous_end
                isempty(current) || append!(current, points[2:end])
                isempty(current) && append!(current, points)
            else
                length(current) >= 3 && push!(components, current)
                current = copy(points)
            end
            previous_end = Int(round(row[8]))
        end
        length(current) >= 3 && push!(components, current)
    end
    # R closes rings explicitly.  Remove a duplicate terminal point before
    # closure so wrapper geometry receives one and only one closing vertex.
    for ring in components
        if length(ring) > 1 && hypot(ring[1][1] - ring[end][1], ring[1][2] - ring[end][2]) <= tol
            pop!(ring)
        end
        push!(ring, ring[1])
    end
    return components
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
