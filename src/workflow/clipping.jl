function _calculate_hull(points::AbstractMatrix, context::RangeContext, spatial::SpatialWorker)
    native = (
        context=spatial.transforms.context,
        forward=spatial.transforms.alpha_forward,
        inverse=spatial.transforms.alpha_inverse,
        geos=spatial.geos,
    )
    return _RB.getDynamicAlphaHull(
        points; fraction=context.fraction, partCount=context.part_count,
        buff=context.buffer_m, initialAlpha=context.initial_alpha,
        coordHeaders=(1, 2), clipToCoast=:no,
        alphaIncrement=context.alpha_increment, alphaCap=context.alpha_cap,
        native,
    )
end

function _cea_projection(hull, spatial::SpatialWorker)
    try
        return as_geos(GO.reproject(
            hull, spatial.transforms.forward; target_crs=EQUAL_AREA_CRS), spatial.geos)
    catch
        return nothing
    end
end

function _interior_terrestrial_projection(hull, spatial::SpatialWorker; projected=nothing)
    started = time_ns()
    try
        projected = isnothing(projected) ? _cea_projection(hull, spatial) : projected
        isnothing(projected) && return nothing
        eligible = try
            !LG.isEmpty(projected, spatial.geos) && LG.isValid(projected, spatial.geos) &&
            LG.covers(spatial.land_prepared, projected, spatial.geos)
        catch
            false
        end
        if eligible
            _s2_stat!(spatial.s2_stats, :interior_covers)
            return projected
        end
        _s2_stat!(spatial.s2_stats, :interior_miss)
        return nothing
    finally
        _s2_time!(spatial.s2_stats, :interior_total, time_ns() - started)
    end
end

function _land_candidate_world(projected, spatial::SpatialWorker)
    spatial.land_candidate_enabled || begin
        _s2_stat!(spatial.s2_stats, :fallback_candidate_disabled)
        return (status=:fallback, candidate=nothing)
    end
    x_span, y_span = try
        (LG.getXMax(projected, spatial.geos) - LG.getXMin(projected, spatial.geos),
         LG.getYMax(projected, spatial.geos) - LG.getYMin(projected, spatial.geos))
    catch
        _s2_stat!(spatial.s2_stats, :fallback_global_extent)
        return (status=:fallback, candidate=nothing)
    end
    (x_span >= REGIONAL_LAND_MAX_FRACTION * _EARTH_WIDTH_M ||
     y_span >= REGIONAL_LAND_MAX_FRACTION * _EARTH_HEIGHT_M) && begin
        _s2_stat!(spatial.s2_stats, :fallback_global_extent)
        return (status=:fallback, candidate=nothing)
    end
    query = try
        started = time_ns()
        value = LG.buffer(LG.envelope(projected, spatial.geos),
                          REGIONAL_LAND_MARGIN_M, 8, spatial.geos)
        _s2_time!(spatial.s2_stats, :regional_buffer, time_ns() - started)
        value
    catch
        _s2_stat!(spatial.s2_stats, :fallback_buffer_failure)
        return (status=:fallback, candidate=nothing)
    end
    candidates = try
        started = time_ns()
        value = LG.query(spatial.land_cea_tree, query; context=spatial.geos)
        _s2_time!(spatial.s2_stats, :regional_tree_query, time_ns() - started)
        value
    catch
        _s2_stat!(spatial.s2_stats, :fallback_query_failure)
        return (status=:fallback, candidate=nothing)
    end
    isempty(candidates) && begin
        _s2_stat!(spatial.s2_stats, :regional_empty_candidates)
        return (status=:empty, candidate=nothing)
    end
    indices = try
        started = time_ns()
        value = sort!(unique!([spatial.land_component_index[c] for c in candidates]))
        _s2_time!(spatial.s2_stats, :regional_index_keys, time_ns() - started)
        value
    catch
        _s2_stat!(spatial.s2_stats, :fallback_index_failure)
        return (status=:fallback, candidate=nothing)
    end
    length(indices) > REGIONAL_LAND_MAX_FRACTION * length(spatial.land_cea_components) && begin
        _s2_stat!(spatial.s2_stats, :fallback_too_many_candidates)
        return (status=:fallback, candidate=nothing)
    end
    started = time_ns()
    key = BitSet(indices)
    hit = haskey(spatial.land_candidate_cache, key)
    _s2_time!(spatial.s2_stats, :regional_cache_lookup, time_ns() - started)
    if hit
        _s2_stat!(spatial.s2_stats, :regional_cache_hit)
        _s2_stat!(spatial.s2_stats, :regional_candidate)
        return (status=:candidate, candidate=spatial.land_candidate_cache[key])
    end
    _s2_stat!(spatial.s2_stats, :regional_cache_miss)
    candidate = let
        started = time_ns()
        selected = spatial.land_wgs84_components[indices]
        world = length(selected) == 1 ? only(selected) :
                GI.Wrappers.GeometryCollection(selected; crs=4326)
        value = (geometry=world, s2=S2.Geography(world))
        _s2_time!(spatial.s2_stats, :regional_candidate_build, time_ns() - started)
        value
    end
    if spatial.candidate_cache_limit > 0
        if length(spatial.land_candidate_cache) >= spatial.candidate_cache_limit
            empty!(spatial.land_candidate_cache)
            _s2_stat!(spatial.s2_stats, :regional_cache_reset)
        end
        spatial.land_candidate_cache[key] = candidate
    end
    _s2_stat!(spatial.s2_stats, :regional_candidate)
    return (status=:candidate, candidate=candidate)
end

function _candidate_terrestrial_projection(hull, projected, spatial::SpatialWorker)
    started = time_ns()
    try
        query_started = time_ns()
        query = _land_candidate_world(projected, spatial)
        _s2_time!(spatial.s2_stats, :regional_index, time_ns() - query_started)
        query.status === :fallback && return nothing
        query.status === :empty && return LG.readgeom("GEOMETRYCOLLECTION EMPTY", spatial.geos)
        candidate = query.candidate
        overlay_started = time_ns()
        clipped = try
            _RB._s2_clip_components(hull, candidate.geometry, :terrestrial;
                                    native=spatial.transforms, world_s2=candidate.s2)
        catch
            _s2_stat!(spatial.s2_stats, :fallback_local_clip_failure)
            return nothing
        finally
            _s2_time!(spatial.s2_stats, :regional_overlay, time_ns() - overlay_started)
        end
        isempty(_RB._polygon_components(clipped)) && begin
            _s2_stat!(spatial.s2_stats, :regional_empty_result)
            return nothing
        end
        try
            valid_geometry(GO.reproject(
                clipped, spatial.transforms.forward; target_crs=EQUAL_AREA_CRS);
                context=spatial.geos)
        catch
            _s2_stat!(spatial.s2_stats, :fallback_local_clip_failure)
            nothing
        end
    finally
        _s2_time!(spatial.s2_stats, :regional_total, time_ns() - started)
    end
end

function clip_hull_to_context(hull, context::RangeContext; spatial=nothing)
    if isnothing(spatial)
        context.clip_to_coast === :no && return valid_geometry(GO.reproject(
            hull; source_crs="EPSG:4326", target_crs=EQUAL_AREA_CRS))
        land_wgs84 = valid_geometry(GO.reproject(
            context.land_cea; source_crs=EQUAL_AREA_CRS, target_crs="EPSG:4326"))
        clipped = _RB._s2_clip_components(hull, land_wgs84, context.clip_to_coast)
        return valid_geometry(GO.reproject(
            clipped; source_crs="EPSG:4326", target_crs=EQUAL_AREA_CRS))
    end
    context.clip_to_coast === :no && return valid_geometry(
        GO.reproject(hull, spatial.transforms.forward; target_crs=EQUAL_AREA_CRS);
        context=spatial.geos)
    if context.clip_to_coast === :terrestrial
        projected = _cea_projection(hull, spatial)
        interior = _interior_terrestrial_projection(hull, spatial; projected)
        isnothing(interior) || return interior
        nonempty = try
            !isnothing(projected) && !LG.isEmpty(projected, spatial.geos)
        catch
            false
        end
        if nonempty
            candidate = _candidate_terrestrial_projection(hull, projected, spatial)
            isnothing(candidate) || return candidate
        end
    end
    _s2_stat!(spatial.s2_stats, :global_overlay)
    started = time_ns()
    try
        clipped = _RB._s2_clip_components(
            hull, spatial.land_wgs84, context.clip_to_coast;
            native=spatial.transforms, world_s2=spatial.land_s2)
        return valid_geometry(GO.reproject(
            clipped, spatial.transforms.forward; target_crs=EQUAL_AREA_CRS);
            context=spatial.geos)
    finally
        _s2_time!(spatial.s2_stats, :global_total, time_ns() - started)
    end
end
