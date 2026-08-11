const _range_land_buffer_cache = Dict{Tuple{Int,Float64},Any}()

function _range_land_coordinates(coords; coordHeaders=(1, 2))
    if coords isa AbstractVector
        length(coords) == 2 || throw(ArgumentError("coordinate vectors must have length two"))
        return reshape(Float64[ismissing(value) ? NaN : value for value in coords], 1, 2)
    end
    return _range_coordinates(coords; coordHeaders)
end

function _proj_crs(crs)
    crs === nothing && throw(ArgumentError("crs must be specified"))
    crs isa Integer && return "EPSG:$(Int(crs))"
    return string(crs)
end

function _buffered_land(scale::Integer, buffer::Real)
    key = (Int(scale), Float64(buffer))
    haskey(_range_land_buffer_cache, key) && return _range_land_buffer_cache[key]
    land = _naturalearth_land(scale)
    isnothing(land) && return nothing
    projected = GeometryOps.reproject(land; source_crs="EPSG:4326", target_crs="+proj=eqearth")
    buffered = GeometryOps.buffer(projected, Float64(buffer))
    _range_land_buffer_cache[key] = buffered
    return buffered
end

"""
    filterByLand(coords; crs=4326, coordHeaders=(1, 2), coastScale=50, buffer=2000)

Return a logical vector indicating whether occurrence coordinates are on land.
The Natural Earth land polygons are buffered by `buffer` metres in Equal Earth
coordinates, matching the two-kilometre coastal tolerance in R `rangeBuilder`.
Missing coordinate rows return `missing`. R uses a bundled GSHHG raster, so
only coastline-adjacent classifications can differ because of source detail.
"""
function filterByLand(coords; crs=4326, coordHeaders=(1, 2),
                      coastScale::Integer=50, buffer::Real=2000)
    coastScale == 50 || throw(ArgumentError("only coastScale=50 is bundled"))
    buffer >= 0 || throw(ArgumentError("buffer must be nonnegative"))
    points = _range_land_coordinates(coords; coordHeaders)
    land = _buffered_land(coastScale, buffer)
    isnothing(land) && return fill(missing, size(points, 1))

    source_crs = _proj_crs(crs)
    wrappers = GeoInterface.Wrappers
    result = Vector{Union{Missing,Bool}}(undef, size(points, 1))
    for i in axes(points, 1)
        if !all(isfinite, view(points, i, :))
            result[i] = missing
            continue
        end
        point = wrappers.Point(points[i, 1], points[i, 2]; crs=crs)
        projected = GeometryOps.reproject(point; source_crs, target_crs="+proj=eqearth")
        result[i] = GeometryOps.intersects(GeometryOps.GEOS(), projected, land)
    end
    return result
end
