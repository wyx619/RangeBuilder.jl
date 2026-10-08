function _build_grid_locator(geometries)
    inverse = Proj.Transformation(EQUAL_AREA_CRS, "EPSG:4326"; always_xy=true)
    cells = Dict{Tuple{Int,Int},Int}()
    try
        for (index, geometry) in enumerate(geometries)
            centroid = LG.centroid(geometry)
            point = GI.Wrappers.Point(
                LG.getGeomX(centroid), LG.getGeomY(centroid); crs=EQUAL_AREA_CRS)
            geographic = GO.reproject(point, inverse; target_crs="EPSG:4326")
            key = _grid_bin(GI.x(geographic), GI.y(geographic))
            isnothing(key) || (cells[key] = index)
        end
    finally
        finalize(inverse)
    end
    return GridLocator(cells, 360, 180)
end

"""Read the package-bundled one-degree map in the CEA working CRS."""
function read_one_degree_map(; target_crs=EQUAL_AREA_CRS)
    target_crs == EQUAL_AREA_CRS ||
        throw(ArgumentError("the bundled workflow grid must use EQUAL_AREA_CRS"))
    payload = _RB._read_one_degree_map_resource()
    ids, geometries = Int[], Any[]
    for (identifier, bytes) in zip(payload["ids"], payload["wkb"])
        raw = LG.readgeom(Vector{Cuchar}(bytes))
        projected = payload["crs"] == target_crs ?
                    valid_geometry(raw) :
                    valid_geometry(GO.reproject(
                        raw; source_crs=payload["crs"], target_crs=target_crs))
        push!(ids, Int(identifier))
        push!(geometries, projected)
    end
    isempty(ids) && error("bundled one-degree map contains no usable cells")
    lookup = IdDict{Any,Int}(geometry => index for (index, geometry) in enumerate(geometries))
    wkb = geometry_wkb.(geometries)
    locator = _build_grid_locator(geometries)
    return GridData(ids, geometries, LG.STRtree(geometries), lookup, wkb, locator)
end

"""Return the package-bundled Natural Earth 1:50m land geometry."""
bundled_land() = _RB._naturalearth_land(50)

function new_range_context(grid::GridData;
                           land=bundled_land(), alpha_min_points=3,
                           buffer_m=10_000, fraction=0.95, part_count=3,
                           initial_alpha=2, alpha_increment=1, alpha_cap=400,
                           clip_to_coast=:terrestrial,
                           candidate_cache_limit=1024)
    alpha_min_points >= 3 || error("alpha_min_points must be at least 3")
    buffer_m > 0 || error("buffer_m must be positive")
    0 < fraction <= 1 || error("fraction must be in (0, 1]")
    part_count >= 1 || error("part_count must be positive")
    initial_alpha >= 0 || error("initial_alpha must be nonnegative")
    alpha_increment > 0 || error("alpha_increment must be positive")
    alpha_cap >= initial_alpha || error("alpha_cap must not be below initial_alpha")
    clip_to_coast in (:no, :terrestrial, :aquatic) || error("invalid clip_to_coast")
    candidate_cache_limit >= 0 || error("candidate_cache_limit must be nonnegative")
    land_cea = valid_geometry(GO.reproject(
        land; source_crs="EPSG:4326", target_crs=EQUAL_AREA_CRS))
    return RangeContext(
        grid, land_cea, Int(alpha_min_points), Float64(buffer_m), Float64(fraction),
        Int(part_count), Float64(initial_alpha), Float64(alpha_increment),
        Float64(alpha_cap), clip_to_coast, Int(candidate_cache_limit),
        geometry_wkb(land_cea), S2ClipStats(),
    )
end

"""Build a context from the bundled grid and bundled land resource."""
prepare_workflow_context(; kwargs...) =
    new_range_context(read_one_degree_map(); kwargs...)

function _workflow_spatial_context()
    context = Proj.proj_context_create()
    forward = Proj.Transformation("EPSG:4326", EQUAL_AREA_CRS;
                                  always_xy=true, ctx=context)
    inverse = Proj.Transformation(EQUAL_AREA_CRS, "EPSG:4326";
                                  always_xy=true, ctx=context)
    alpha_forward = Proj.Transformation("EPSG:4326", "+proj=eqearth";
                                         always_xy=true, ctx=context)
    alpha_inverse = Proj.Transformation("+proj=eqearth", "EPSG:4326";
                                         always_xy=true, ctx=context)
    holder = WorkflowSpatialContext(
        context, forward, inverse, alpha_forward, alpha_inverse, LG.GEOSContext())
    finalizer(holder) do value
        finalize(value.forward)
        finalize(value.inverse)
        finalize(value.alpha_forward)
        finalize(value.alpha_inverse)
        finalize(value.geos)
        value.context != C_NULL && Proj.proj_context_destroy(value.context)
        value.context = C_NULL
    end
    return holder
end
