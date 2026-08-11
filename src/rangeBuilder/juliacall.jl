"""Serialize a GeoInterface polygon or multipolygon as WKT.

This small text-only adapter is intended for foreign-language callers.  It
keeps GeoInterface wrapper objects inside Julia, where the geometry kernel
operates, and returns a format that R `terra`/`sf` can read directly.
"""
function _geometry_wkt(geometry)
    trait = GeoInterface.geomtrait(geometry)
    coordinates = GeoInterface.coordinates(trait, geometry)
    position(point) = string(Float64(point[1]), " ", Float64(point[2]))
    ring_wkt(ring) = string("(", join(position.(ring), ", "), ")")
    polygon_wkt(rings) = string("(", join(ring_wkt.(rings), ", "), ")")

    if trait isa GeoInterface.PolygonTrait
        return string("POLYGON ", polygon_wkt(coordinates))
    elseif trait isa GeoInterface.MultiPolygonTrait
        return string("MULTIPOLYGON (", join(polygon_wkt.(coordinates), ", "), ")")
    end
    throw(ArgumentError("expected a polygon or multipolygon, got $(typeof(geometry))"))
end

"""
    getDynamicAlphaHullWKT(x; kwargs...)

Run `getDynamicAlphaHull` and return the resulting polygon as WKT.  This is
the text bridge used by R via JuliaCall; all alpha-hull computation remains in
Julia.  Keyword arguments are exactly those of `getDynamicAlphaHull`.
"""
function getDynamicAlphaHullWKT(x; kwargs...)
    result = getDynamicAlphaHull(x; kwargs...)
    return _geometry_wkt(result.hull)
end
