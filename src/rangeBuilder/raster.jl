"""Return the combined GeoInterface extent of a collection of geometries."""
function getExtentOfList(shapes)
    isempty(shapes) && throw(ArgumentError("shapes must not be empty"))
    extents = [GeoInterface.extent(shape) for shape in shapes if !isnothing(shape)]
    isempty(extents) && throw(ArgumentError("shapes contains no non-empty geometries"))
    xmin = minimum(ext.X[1] for ext in extents)
    xmax = maximum(ext.X[2] for ext in extents)
    ymin = minimum(ext.Y[1] for ext in extents)
    ymax = maximum(ext.Y[2] for ext in extents)
    return Rasters.Extent(X=(Float64(xmin), Float64(xmax)),
                          Y=(Float64(ymin), Float64(ymax)))
end

function _named_polygon_pairs(polyList, speciesNames)
    if polyList isa NamedTuple
        return collect(pairs(polyList))
    elseif polyList isa AbstractDict
        isempty(polyList) && throw(ArgumentError("polyList must not be empty"))
        return Pair{Symbol,Any}[Symbol(k) => v for (k, v) in polyList]
    elseif polyList isa AbstractVector
        isnothing(speciesNames) &&
            throw(ArgumentError("polyList must be named; pass speciesNames for a vector"))
        length(speciesNames) == length(polyList) ||
            throw(ArgumentError("speciesNames and polyList must have equal length"))
        return Pair{Symbol,Any}[Symbol(k) => v for (k, v) in zip(speciesNames, polyList)]
    end
    throw(ArgumentError("polyList must be a NamedTuple, dictionary, or vector"))
end

function _raster_crs(crs)
    crs === nothing && return nothing
    crs isa Integer && return Rasters.EPSG(Int(crs))
    crs isa AbstractString && return try
        Rasters.ProjString(String(crs))
    catch
        nothing
    end
    crs isa Rasters.EPSG || crs isa Rasters.ProjString || return nothing
    return crs
end

function _is_longlat_crs(crs)
    crs == 4326 || crs == Rasters.EPSG(4326) ||
        (crs isa AbstractString && (crs == "EPSG:4326" || occursin("longlat", lowercase(crs))))
end

function _geometry_coordinates!(out, value)
    if value isa Tuple || value isa AbstractVector
        if length(value) >= 2 && value[1] isa Real && value[2] isa Real
            push!(out, (Float64(value[1]), Float64(value[2])))
        else
            for child in value
                _geometry_coordinates!(out, child)
            end
        end
    end
    return out
end

function _fallback_interior_point(poly)
    try
        candidate = GeometryOps.centroid(poly)
        candidate isa Tuple && length(candidate) >= 2 &&
            return (Float64(candidate[1]), Float64(candidate[2]))
    catch
    end
    coords = Tuple{Float64,Float64}[]
    _geometry_coordinates!(coords, GeoInterface.coordinates(poly))
    isempty(coords) && return nothing
    return (sum(first, coords) / length(coords), sum(last, coords) / length(coords))
end

function _retain_small_range!(raster, poly, grid_extent, resolution)
    point = _fallback_interior_point(poly)
    isnothing(point) && return raster
    xmin, ymin = grid_extent.X[1], grid_extent.Y[1]
    ix = clamp(floor(Int, (point[1] - xmin) / resolution) + 1, 1, size(raster, 1))
    iy = clamp(floor(Int, (point[2] - ymin) / resolution) + 1, 1, size(raster, 2))
    raster[ix, iy] = 1
    return raster
end

"""
    rasterStackFromPolyList(polyList; resolution=50000, retainSmallRanges=true,
                            extent=:auto, speciesNames=nothing)

Create a `Rasters.RasterStack` with one `missing`/`1` layer per named polygon.
Polygon boundaries use cell centers, matching `terra::rasterize` in the R
`rangeBuilder` package. `polyList` may be a `NamedTuple`, dictionary, or vector
with `speciesNames` supplied.
"""
function rasterStackFromPolyList(polyList; resolution=50000,
                                 retainSmallRanges::Bool=true, extent=:auto,
                                 speciesNames=nothing, threaded::Bool=true)
    resolution isa Real && resolution > 0 ||
        throw(ArgumentError("resolution must be a positive number"))
    pairs_ = _named_polygon_pairs(polyList, speciesNames)
    names_ = first.(pairs_)
    length(unique(names_)) == length(names_) || throw(ArgumentError("species names must be unique"))

    nonempty = Pair{Symbol,Any}[p for p in pairs_ if !isnothing(p.second)]
    isempty(nonempty) && throw(ArgumentError("polyList contains no non-empty geometries"))
    first_poly = first(nonempty).second
    crs_value = GeoInterface.crs(first_poly)
    longlat = _is_longlat_crs(crs_value)
    longlat && resolution > 90 &&
        throw(ArgumentError("input data are in longlat, therefore resolution must be in decimal degrees"))
    !longlat && crs_value !== nothing && resolution < 90 &&
        throw(ArgumentError("input data are projected, but resolution is unexpectedly small"))

    for p in nonempty
        GeoInterface.crs(p.second) == crs_value ||
            throw(ArgumentError("all polygons must use the same CRS"))
    end
    master_extent = if extent === :auto || extent == "auto"
        getExtentOfList(last.(nonempty))
    elseif extent isa AbstractVector && length(extent) == 4 && all(x -> x isa Real, extent)
        Rasters.Extent(X=(Float64(extent[1]), Float64(extent[2])),
                       Y=(Float64(extent[3]), Float64(extent[4])))
    else
        throw(ArgumentError("extent must be :auto or [minLong, maxLong, minLat, maxLat]"))
    end

    raster_crs = _raster_crs(crs_value)
    # Terra expands the upper edge when the requested extent is not an exact
    # multiple of the resolution. Build the same regular cell-center grid.
    nx = max(1, ceil(Int, (master_extent.X[2] - master_extent.X[1]) / resolution))
    ny = max(1, ceil(Int, (master_extent.Y[2] - master_extent.Y[1]) / resolution))
    grid_extent = Rasters.Extent(
        X=(master_extent.X[1], master_extent.X[1] + nx * Float64(resolution)),
        Y=(master_extent.Y[1], master_extent.Y[1] + ny * Float64(resolution)),
    )
    layers = Pair{Symbol,Any}[]
    for (name, poly) in pairs_
        isnothing(poly) && continue
        template = Rasters.Raster(grid_extent; size=(nx, ny),
                                  sampling=Rasters.Intervals(Rasters.Center()),
                                  crs=raster_crs, missingval=missing)
        layer = Rasters.rasterize(last, poly; to=template, shape=:polygon,
                                  boundary=:center, fill=1, missingval=missing,
                                  threaded)
        if retainSmallRanges && all(ismissing, parent(layer))
            _retain_small_range!(layer, poly, grid_extent, Float64(resolution))
        elseif !retainSmallRanges && all(ismissing, parent(layer))
            continue
        end
        push!(layers, name => layer)
    end
    isempty(layers) && throw(ArgumentError("no polygon intersects the raster extent"))
    keys_ = Tuple(first.(layers))
    vals_ = Tuple(last.(layers))
    return Rasters.RasterStack(NamedTuple{keys_}(vals_))
end

"""
    speciesRichness(stack; zeroToMissing=true, name=:richness)

Count the present species layers in a `Rasters.RasterStack`. This mirrors the
R workflow `terra::app(stack, sum, na.rm=TRUE)` followed by converting zeroes
to `NA`. Layers are expected to contain `1` for presence and `missing` for
absence, as produced by [`rasterStackFromPolyList`](@ref).
"""
function speciesRichness(stack::Rasters.AbstractRasterStack;
                         zeroToMissing::Bool=true, name::Symbol=:richness)
    layer_names = keys(stack)
    isempty(layer_names) && throw(ArgumentError("stack must contain at least one layer"))
    first_layer = stack[first(layer_names)]
    counts = zeros(Int, size(first_layer))
    present = falses(size(first_layer))
    for layer_name in layer_names
        layer = stack[layer_name]
        size(layer) == size(first_layer) ||
            throw(ArgumentError("all raster layers must have equal dimensions"))
        for index in eachindex(counts)
            value = layer[index]
            if !ismissing(value)
                counts[index] += 1
                present[index] = true
            end
        end
    end
    if zeroToMissing
        data = Array{Union{Missing,Int}}(undef, size(counts))
        for index in eachindex(data)
            data[index] = present[index] ? counts[index] : missing
        end
    else
        data = counts
    end
    return Rasters.Raster(data, Rasters.dims(first_layer); missingval=missing, name)
end
