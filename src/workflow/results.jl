function _result(task::SpeciesTask, method; grid_ids=Int[], area_km2=missing,
                 area_status=:failed, geometry=nothing, error_message=nothing,
                 status=nothing)
    success = isnothing(error_message) && !isempty(grid_ids) &&
              !ismissing(area_km2) && isfinite(area_km2) && area_km2 > 0
    final_status = isnothing(status) ? (success ? :success : :failed) : status
    return SpeciesRangeResult(
        task.name, size(task.points, 1), method, final_status, success,
        task.metadata, task.source_rows, Int.(grid_ids),
        success ? Float64(area_km2) : missing, area_status,
        success ? geometry : nothing,
        success ? nothing : something(error_message, "no valid grid or positive area"),
    )
end

function _alpha_result(task::SpeciesTask, hull_value, context, spatial)
    method = Symbol("alpha_hull_$(context.buffer_m / 1000)km_buffer_clipped_to_land")
    hull_value isa AbstractString &&
        return _result(task, method; error_message=hull_value)
    try
        geometry = clip_hull_to_context(hull_value.hull, context; spatial)
        geometry = valid_geometry(geometry; context=spatial.geos)
        ids = intersecting_grid_ids(geometry, context.grid; spatial)
        return _result(task, method; grid_ids=ids,
                       area_km2=_ellipsoidal_area_km2(geometry, spatial),
                       area_status=:ok, geometry)
    catch error
        return _result(task, method; error_message=sprint(showerror, error))
    end
end

function _direct_result(task::SpeciesTask, context, spatial;
                        method=Symbol("direct_point_grid_plus_$(context.buffer_m / 1000)km_buffer_area"))
    try
        ids = _point_grid_ids(task.points, context.grid, spatial)
        area = _point_buffer_area(task.points, context, spatial)
        return _result(task, method; grid_ids=ids, area_km2=area.area_km2,
                       area_status=area.status, geometry=area.geometry)
    catch error
        return _result(task, method; error_message=sprint(showerror, error))
    end
end

function _excluded_result(task::SpeciesTask, status::Symbol)
    method = Symbol("excluded_$(status)")
    return _result(task, method; status)
end

function _metadata_value(values, rows)
    selected = Any[values[row] for row in rows if !ismissing(values[row])]
    isempty(selected) && return missing
    unique_values = unique(selected)
    return length(unique_values) == 1 ? only(unique_values) : unique_values
end

function _occurrence_tasks(input::OccurrenceInput)
    groups = Dict{Symbol,Vector{Int}}()
    for (row, species) in enumerate(input.species)
        push!(get!(groups, Symbol(species), Int[]), row)
    end
    tasks = SpeciesTask[]
    metadata_names = propertynames(input.metadata)
    for name in sort!(collect(keys(groups)))
        rows = groups[name]
        points = hcat(input.longitude[rows], input.latitude[rows])
        metadata = NamedTuple{metadata_names}(
            Tuple(_metadata_value(getproperty(input.metadata, field), rows)
                  for field in metadata_names))
        push!(tasks, SpeciesTask(
            name, Matrix{Float64}(points), metadata, input.source_rows[rows]))
    end
    return tasks
end

function _matrix_tasks(species_points)
    tasks = SpeciesTask[]
    for (name, points) in _species_pairs(species_points)
        cleaned = clean_species_points(points)
        push!(tasks, SpeciesTask(name, Matrix{Float64}(cleaned), NamedTuple(), Int[]))
    end
    sort!(tasks; by=task -> task.name)
    return tasks
end

function _build_species_ranges(tasks::Vector{SpeciesTask}, context::RangeContext;
                                workers::Integer)
    started = time_ns()
    direct = [task for task in tasks if size(task.points, 1) < context.alpha_min_points]
    # Keep all sufficiently large point sets on the R-compatible SHull path.
    # Collinear triples can still recover through the same alpha/MCH fallback
    # as the reference workflow, so a pre-emptive geometric exclusion would
    # change accepted results.
    alpha = [task for task in tasks if size(task.points, 1) >= context.alpha_min_points]
    results = SpeciesRangeResult[]
    append!(results, _map_direct_ranges(direct, context; workers))
    append!(results, _map_ranges(alpha, context; workers))
    sort!(results; by=result -> result.species)
    elapsed = (time_ns() - started) / 1e9
    return SpeciesRangesResult(results, elapsed,
                               s2_clip_stats(context), s2_clip_timings(context))
end

"""
    build_species_ranges(input::OccurrenceInput, context; workers=Threads.nthreads())

Build ranges from standardized table columns while carrying species metadata
and source-row indices into each `SpeciesRangeResult`.
"""
build_species_ranges(input::OccurrenceInput, context::RangeContext;
                     workers::Integer=nthreads()) =
    _build_species_ranges(_occurrence_tasks(input), context; workers)

"""
    build_species_ranges(species_points, context; workers=Threads.nthreads())

Compatibility entry point for a NamedTuple or dictionary of species matrices.
"""
build_species_ranges(species_points, context::RangeContext;
                     workers::Integer=nthreads()) =
    _build_species_ranges(_matrix_tasks(species_points), context; workers)
