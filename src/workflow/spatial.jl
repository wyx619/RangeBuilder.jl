function _spatial_worker(context::RangeContext)
    tls = task_local_storage()
    key = (:range_builder_workflow_spatial_worker, objectid(context))
    worker = get(tls, key, nothing)
    worker === nothing || return worker
    transforms = _workflow_spatial_context()
    geos = transforms.geos
    land = LG.readgeom(context.land_wkb, geos)
    land_wgs84 = valid_geometry(
        GO.reproject(land, transforms.inverse; target_crs="EPSG:4326"); context=geos)
    land_s2 = S2.Geography(land_wgs84)
    land_prepared = LG.prepareGeom(land, geos)
    land_cea_components = _RB._polygon_components(land)
    land_wgs84_components = _RB._polygon_components(land_wgs84)
    enabled = length(land_cea_components) == length(land_wgs84_components) &&
              !isempty(land_cea_components)
    land_tree = LG.STRtree(land_cea_components; context=geos)
    component_index = IdDict{Any,Int}(
        geometry => index for (index, geometry) in enumerate(land_cea_components))
    grid_geometries = [LG.readgeom(wkb, geos) for wkb in context.grid.wkb]
    grid_tree = LG.STRtree(grid_geometries; context=geos)
    grid_index = IdDict{Any,Int}(
        geometry => index for (index, geometry) in enumerate(grid_geometries))
    worker = SpatialWorker(
        transforms, geos, land, land_wgs84, land_s2, land_prepared,
        land_cea_components, land_wgs84_components, land_tree, component_index,
        Dict{BitSet,Any}(), context.candidate_cache_limit, enabled,
        context.s2_stats, grid_geometries, grid_tree, grid_index,
    )
    tls[key] = worker
    return worker
end

function _direct_worker(context::RangeContext)
    tls = task_local_storage()
    key = (:range_builder_workflow_direct_worker, objectid(context))
    worker = get(tls, key, nothing)
    worker === nothing || return worker
    transforms = _workflow_spatial_context()
    geos = transforms.geos
    land = LG.readgeom(context.land_wkb, geos)
    worker = DirectWorker(
        transforms, geos, land, LG.prepareGeom(land, geos),
        nothing, nothing, Dict{Int,Any}())
    tls[key] = worker
    return worker
end

function _locator_candidates(locator::GridLocator, key::Tuple{Int,Int})
    lon, lat = key
    values = Int[]
    for dlat in -1:1, dlon in -1:1
        candidate = (lon + dlon, lat + dlat)
        1 <= candidate[1] <= locator.lon_bins && 1 <= candidate[2] <= locator.lat_bins || continue
        index = get(locator.cells, candidate, nothing)
        isnothing(index) || push!(values, index)
    end
    return unique!(values)
end

function _direct_candidate_geometries(spatial::DirectWorker, grid::GridData, indices)
    for index in indices
        haskey(spatial.candidate_geometries, index) ||
            (spatial.candidate_geometries[index] = LG.readgeom(grid.wkb[index], spatial.geos))
    end
    return spatial.candidate_geometries
end

function _direct_grid_data(spatial::DirectWorker, grid::GridData)
    isnothing(spatial.grid_geometries) || return spatial.grid_geometries, spatial.grid_tree
    geometries = [LG.readgeom(wkb, spatial.geos) for wkb in grid.wkb]
    spatial.grid_geometries = geometries
    spatial.grid_tree = LG.STRtree(geometries; context=spatial.geos)
    return geometries, spatial.grid_tree
end

function nearest_grid_index(point, geometries, geos::LG.GEOSContext)
    distances = [LG.distance(point, candidate, geos) for candidate in geometries]
    return argmin(distances)
end

function _point_grid_ids(points::AbstractMatrix, grid::GridData, spatial::DirectWorker)
    identifiers = Int[]
    for row in eachrow(points)
        longitude, latitude = Float64(row[1]), Float64(row[2])
        key = _grid_bin(longitude, latitude)
        point_wgs84 = GI.Wrappers.Point(longitude, latitude; crs=4326)
        point = as_geos(GO.reproject(
            point_wgs84, spatial.transforms.forward; target_crs=EQUAL_AREA_CRS), spatial.geos)
        candidate_indices = isnothing(key) || isnothing(grid.locator) ? Int[] :
            _locator_candidates(grid.locator, key)
        if !isempty(candidate_indices)
            geometries = _direct_candidate_geometries(spatial, grid, candidate_indices)
            contained = findfirst(
                index -> LG.within(point, geometries[index], spatial.geos), candidate_indices)
            isnothing(contained) || begin
                push!(identifiers, grid.ids[candidate_indices[contained]])
                continue
            end
        end
        geometries, tree = _direct_grid_data(spatial, grid)
        candidates = LG.query(tree, point; context=spatial.geos)
        contained = findfirst(
            geometry -> LG.within(point, geometry, spatial.geos), candidates)
        if !isnothing(contained)
            push!(identifiers, grid.ids[findfirst(==(candidates[contained]), geometries)])
        else
            push!(identifiers, grid.ids[nearest_grid_index(point, geometries, spatial.geos)])
        end
    end
    return sort!(unique!(identifiers))
end

function intersecting_grid_ids(geometry, grid::GridData; spatial=nothing)
    if isnothing(spatial)
        polygon = valid_geometry(geometry)
        tree, lookup, geos = grid.tree, grid.geometry_index, nothing
    else
        geos = spatial.geos
        polygon = valid_geometry(geometry; context=geos)
        tree, lookup = spatial.tree, spatial.geometry_index
    end
    context = isnothing(geos) ? LG.get_context(polygon) : geos
    candidates = LG.query(tree, polygon; context)
    ids = Int[]
    for candidate in candidates
        LG.intersects(candidate, polygon, context) && push!(ids, grid.ids[lookup[candidate]])
    end
    return sort!(unique!(ids))
end

function _point_buffer_area(points::AbstractMatrix, context::RangeContext, spatial::DirectWorker)
    latitude = sum(points[:, 2]) / size(points, 1)
    longitude = sum(points[:, 1]) / size(points, 1)
    local_laea = "+proj=laea +lat_0=$(latitude) +lon_0=$(longitude) +datum=WGS84 +units=m +no_defs"
    to_laea = Proj.Transformation("EPSG:4326", local_laea;
                                  always_xy=true, ctx=spatial.transforms.context)
    to_cea = Proj.Transformation(local_laea, EQUAL_AREA_CRS;
                                 always_xy=true, ctx=spatial.transforms.context)
    try
        buffers = LG.Geometry[]
        for row in eachrow(points)
            point = GI.Wrappers.Point(Float64(row[1]), Float64(row[2]); crs=4326)
            point_laea = GO.reproject(point, to_laea; target_crs=local_laea)
            push!(buffers, LG.buffer(as_geos(point_laea, spatial.geos),
                                     context.buffer_m, BUFFER_QUADRANT_SEGMENTS, spatial.geos))
        end
        buffered_laea = reduce((left, right) -> LG.union(left, right, spatial.geos), buffers)
        buffered = as_geos(GO.reproject(
            buffered_laea, to_cea; target_crs=EQUAL_AREA_CRS), spatial.geos)
        clipped = try
            LG.covers(spatial.land_prepared, buffered, spatial.geos) ? buffered :
                LG.intersection(buffered, spatial.land_cea, spatial.geos)
        catch
            nothing
        end
        valid_clip = !isnothing(clipped) && LG.area(clipped, spatial.geos) > 0
        geometry = valid_clip ? clipped : buffered
        return (
            geometry=geometry,
            area_km2=_ellipsoidal_area_km2(geometry, spatial),
            status=valid_clip ? :ok : :ok_unclipped_fallback,
        )
    finally
        finalize(to_laea)
        finalize(to_cea)
    end
end

function _ellipsoidal_area_km2(geometry, spatial::SpatialWorker)
    geographic = GO.reproject(geometry, spatial.transforms.inverse; target_crs="EPSG:4326")
    return abs(GO.area(GO.Geodesic(), geographic)) / 1e6
end

function _ellipsoidal_area_km2(geometry, spatial::DirectWorker)
    geographic = GO.reproject(geometry, spatial.transforms.inverse; target_crs="EPSG:4326")
    return abs(GO.area(GO.Geodesic(), geographic)) / 1e6
end
