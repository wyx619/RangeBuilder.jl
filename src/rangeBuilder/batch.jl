function _named_species_pairs(speciesPoints)
    if speciesPoints isa NamedTuple
        return Pair{Symbol,Any}[Symbol(name) => points for (name, points) in pairs(speciesPoints)]
    elseif speciesPoints isa AbstractDict
        return Pair{Symbol,Any}[Symbol(name) => points for (name, points) in speciesPoints]
    end
    throw(ArgumentError("speciesPoints must be a NamedTuple or dictionary keyed by species name"))
end

function _has_2d_support(points::AbstractMatrix)
    size(points, 1) >= 3 || return false
    origin = view(points, 1, :)
    scale = max(1.0, maximum(abs, points))
    tolerance = 64eps(Float64) * scale^2
    for i in 2:size(points, 1)-1
        dx, dy = points[i, 1] - origin[1], points[i, 2] - origin[2]
        dx == 0 && dy == 0 && continue
        for j in i + 1:size(points, 1)
            cross = dx * (points[j, 2] - origin[2]) - dy * (points[j, 1] - origin[1])
            abs(cross) > tolerance && return true
        end
    end
    return false
end

"""
    buildRanges(speciesPoints; coordHeaders=(1, 2), kwargs...)

Build conservative alpha-hull ranges for a named collection of species point
sets. The return value is `(ranges, excluded)`: `ranges` maps successfully
modelled species to GeoInterface polygons, while `excluded` records a status
and unique valid point count for every unmodelled species. One- and two-point
sets are never converted into artificial buffered ranges.
"""
function buildRanges(speciesPoints; coordHeaders=(1, 2), kwargs...)
    ranges = Dict{Symbol,Any}()
    excluded = Dict{Symbol,NamedTuple{(:status, :npoints),Tuple{Symbol,Int}}}()
    for (name, observations) in _named_species_pairs(speciesPoints)
        points = try
            _unique_finite_coordinates(_range_coordinates(observations; coordHeaders))
        catch
            zeros(Float64, 0, 2)
        end
        npoints = size(points, 1)
        if npoints < 3
            excluded[name] = (status=:insufficient_points, npoints=npoints)
        elseif !_has_2d_support(points)
            excluded[name] = (status=:degenerate_geometry, npoints=npoints)
        else
            result = try
                getDynamicAlphaHull(observations; coordHeaders, kwargs...)
            catch
                nothing
            end
            if isnothing(result)
                excluded[name] = (status=:construction_failed, npoints=npoints)
            else
                ranges[name] = result.hull
            end
        end
    end
    return (ranges=ranges, excluded=excluded)
end
