function _normalize_species(value)
    value === missing && return ""
    return strip(replace(replace(string(value), "_" => " "), r"\s+" => " "))
end

function _as_float(value)
    value === missing && return NaN
    try
        return Float64(value)
    catch
        try
            return parse(Float64, string(value))
        catch
            return NaN
        end
    end
end

function _table_columns(table)
    columns = Tables.columns(table)
    names = Symbol.(Tables.columnnames(columns))
    return columns, names
end

function _table_column(table, selector, nrows; required=true, label="column")
    selector === nothing && begin
        required && throw(ArgumentError("$label must be specified"))
        return fill(missing, nrows)
    end
    if selector isa AbstractVector
        values = collect(selector)
    else
        columns, names = _table_columns(table)
        name = selector isa Symbol ? selector : Symbol(selector)
        name in names || throw(ArgumentError("$label $name is not present in the table"))
        values = collect(Tables.getcolumn(columns, name))
    end
    length(values) == nrows ||
        throw(ArgumentError("$label has $(length(values)) rows; expected $nrows"))
    return values
end

function _metadata_mappings(table, metadata, kwargs, core_names)
    mappings = Pair{Symbol,Any}[]
    if metadata === :all
        _, names = _table_columns(table)
        for name in names
            name in core_names || push!(mappings, name => name)
        end
    elseif metadata isa NamedTuple
        append!(mappings, [Symbol(name) => value for (name, value) in pairs(metadata)])
    elseif metadata isa AbstractDict
        append!(mappings, [Symbol(name) => value for (name, value) in metadata])
    elseif metadata isa Tuple || metadata isa AbstractVector
        append!(mappings, [Symbol(name) => name for name in metadata])
    elseif metadata !== nothing
        throw(ArgumentError("metadata must be :all, a NamedTuple, dictionary, tuple, or nothing"))
    end
    for (name, selector) in pairs(kwargs)
        filter!(mapping -> first(mapping) != name, mappings)
        push!(mappings, name => selector)
    end
    names = first.(mappings)
    length(unique(names)) == length(names) || throw(ArgumentError("metadata names must be unique"))
    any(name -> name in core_names, names) &&
        throw(ArgumentError("metadata names cannot overwrite species/longitude/latitude"))
    return mappings
end

"""
    standardize_occurrences(table; species=:species, longitude=:longitude,
                            latitude=:latitude, coordinate_digits=2,
                            metadata=NamedTuple(), kwargs...)

Convert any Tables.jl-compatible table into the workflow's standard occurrence
input. Column selectors may be symbols or strings; `nothing` creates an all-
missing optional column. Coordinates are rounded to `coordinate_digits` decimal
places before de-duplication (default `2`); the retained coordinates are the
rounded values. Additional keyword arguments are metadata mappings, so
`family=:Accepted_family, order=nothing` is supported directly.
"""
function standardize_occurrences(table; species=:species, longitude=:longitude,
                                 latitude=:latitude, coordinate_digits::Integer=2,
                                 metadata=NamedTuple(), kwargs...)
    coordinate_digits >= 0 ||
        throw(ArgumentError("coordinate_digits must be nonnegative"))
    columns, names = _table_columns(table)
    species_name = species isa Symbol ? species : Symbol(species)
    species_name in names || throw(ArgumentError("species column $species_name is not present"))
    raw_species = collect(Tables.getcolumn(columns, species_name))
    nrows = length(raw_species)
    raw_longitude = _table_column(table, longitude, nrows; label="longitude column")
    raw_latitude = _table_column(table, latitude, nrows; label="latitude column")
    core_names = Set((species_name,
                      longitude isa Symbol ? longitude : Symbol(longitude),
                      latitude isa Symbol ? latitude : Symbol(latitude)))
    mappings = _metadata_mappings(table, metadata, kwargs, core_names)
    raw_metadata = [(name, _table_column(table, selector, nrows;
                                         required=false, label="metadata column"))
                    for (name, selector) in mappings]

    species_values = String[]
    longitude_values = Float64[]
    latitude_values = Float64[]
    source_rows = Int[]
    metadata_values = [Any[] for _ in raw_metadata]
    seen = Set{Tuple{String,Float64,Float64}}()
    for row in 1:nrows
        species_value = _normalize_species(raw_species[row])
        longitude_value = _as_float(raw_longitude[row])
        latitude_value = _as_float(raw_latitude[row])
        isfinite(longitude_value) && isfinite(latitude_value) && !isempty(species_value) || continue
        -180 <= longitude_value <= 180 && -90 <= latitude_value <= 90 || continue
        longitude_value = round(longitude_value; digits=coordinate_digits)
        latitude_value = round(latitude_value; digits=coordinate_digits)
        key = (species_value, longitude_value, latitude_value)
        key in seen && continue
        push!(seen, key)
        push!(species_values, species_value)
        push!(longitude_values, longitude_value)
        push!(latitude_values, latitude_value)
        push!(source_rows, row)
        for (index, (_, values)) in enumerate(raw_metadata)
            push!(metadata_values[index], values[row])
        end
    end
    metadata_names = Tuple(first.(raw_metadata))
    metadata_tuple = NamedTuple{metadata_names}(Tuple(metadata_values))
    return OccurrenceInput(
        species_values, longitude_values, latitude_values, metadata_tuple, source_rows)
end

function clean_species_points(points)
    matrix = points isa AbstractMatrix ? points : Matrix{Float64}(points)
    size(matrix, 2) == 2 || throw(ArgumentError("species points must be an n×2 matrix"))
    finite = [all(isfinite, matrix[i, :]) &&
              -180 <= matrix[i, 1] <= 180 && -90 <= matrix[i, 2] <= 90
              for i in axes(matrix, 1)]
    return unique(matrix[finite, :]; dims=1)
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

@inline function _shull_precision_degenerate(error)
    error isa ArgumentError || return false
    return occursin("at least three points remain after SHull duplicate removal",
                    sprint(showerror, error))
end

function _species_pairs(species_points)
    if species_points isa NamedTuple
        return [(Symbol(name), value) for (name, value) in pairs(species_points)]
    elseif species_points isa AbstractDict
        return [(Symbol(name), value) for (name, value) in species_points]
    end
    throw(ArgumentError("species_points must be a NamedTuple or dictionary"))
end

function _map_ranges(tasks, context::RangeContext; workers::Integer)
    results = Vector{Union{Nothing,SpeciesRangeResult}}(nothing, length(tasks))
    isempty(tasks) && return SpeciesRangeResult[]
    1 <= workers <= nthreads() ||
        throw(ArgumentError("workers must be between 1 and $(nthreads())"))
    next_index = Atomic{Int}(1)
    @sync for _ in 1:min(Int(workers), length(tasks))
        Threads.@spawn begin
            spatial = _spatial_worker(context)
            direct = nothing
            while true
                index = atomic_add!(next_index, 1)
                index > length(tasks) && break
                task = tasks[index]
                hull_value = try
                    _calculate_hull(task.points, context, spatial)
                catch error
                    _shull_precision_degenerate(error) ? :shull_precision_degenerate :
                                                        sprint(showerror, error)
                end
                if hull_value === :shull_precision_degenerate
                    # SHull uses Float32 topology, so distinct Float64 input
                    # coordinates can collapse to fewer than three vertices.
                    # Do not manufacture an alpha hull in this case: retain
                    # the original points and use the established direct
                    # point-grid plus buffer-area semantics instead.
                    isnothing(direct) && (direct = _direct_worker(context))
                    results[index] = _direct_result(
                        task, context, direct;
                        method=Symbol(
                            "direct_point_grid_plus_$(context.buffer_m / 1000)km_" *
                            "buffer_area_shull_precision_fallback"),
                    )
                else
                    results[index] = _alpha_result(task, hull_value, context, spatial)
                end
            end
        end
    end
    return SpeciesRangeResult[r for r in results]
end

function _map_direct_ranges(tasks, context::RangeContext; workers::Integer)
    results = Vector{Union{Nothing,SpeciesRangeResult}}(nothing, length(tasks))
    isempty(tasks) && return SpeciesRangeResult[]
    1 <= workers <= nthreads() ||
        throw(ArgumentError("workers must be between 1 and $(nthreads())"))
    next_index = Atomic{Int}(1)
    @sync for _ in 1:min(Int(workers), length(tasks))
        Threads.@spawn begin
            spatial = _direct_worker(context)
            while true
                index = atomic_add!(next_index, 1)
                index > length(tasks) && break
                task = tasks[index]
                results[index] = _direct_result(task, context, spatial)
            end
        end
    end
    return SpeciesRangeResult[r for r in results]
end
