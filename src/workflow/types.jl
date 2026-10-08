as_geos(geometry) = geometry isa LG.Geometry ? geometry : GI.convert(LG, geometry)

function as_geos(geometry, context::LG.GEOSContext)
    geometry isa LG.Geometry && return geometry
    trait = GI.geomtrait(geometry)
    return GI.convert(LG.geointerface_geomtype(trait), trait, geometry; context)
end

function valid_geometry(geometry; context=nothing)
    value = isnothing(context) ? as_geos(geometry) : as_geos(geometry, context)
    geos = isnothing(context) ? LG.get_context(value) : context
    try
        return LG.isValid(value, geos) ? value : LG.makeValid(value, geos)
    catch
        return LG.makeValid(value, geos)
    end
end

function geometry_wkb(geometry)
    value = as_geos(geometry)
    geos = LG.get_context(value)
    return Vector{Cuchar}(LG.writegeom(value, LG.WKBWriter(geos), geos))
end

function union_all(geometries; context=nothing)
    isempty(geometries) && return nothing
    if isnothing(context)
        return reduce((left, right) -> GO.union(GO.GEOS(), left, right), geometries)
    end
    return reduce((left, right) -> LG.union(left, right, context), geometries)
end

struct GridData
    ids::Vector{Int}
    geometries::Vector{Any}
    tree::Any
    geometry_index::IdDict{Any,Int}
    wkb::Vector{Vector{Cuchar}}
    locator::Any
end

GridData(ids, geometries, tree, geometry_index) =
    GridData(ids, geometries, tree, geometry_index, Vector{Vector{Cuchar}}(), nothing)

struct GridLocator
    cells::Dict{Tuple{Int,Int},Int}
    lon_bins::Int
    lat_bins::Int
end

@inline function _grid_bin(longitude::Real, latitude::Real)
    isfinite(longitude) && isfinite(latitude) || return nothing
    (-180 <= longitude <= 180 && -90 <= latitude <= 90) || return nothing
    lon = longitude == 180 ? 360 : floor(Int, longitude + 180) + 1
    lat = latitude == 90 ? 180 : floor(Int, latitude + 90) + 1
    return (clamp(lon, 1, 360), clamp(lat, 1, 180))
end

struct RangeContext
    grid::GridData
    land_cea
    alpha_min_points::Int
    buffer_m::Float64
    fraction::Float64
    part_count::Int
    initial_alpha::Float64
    alpha_increment::Float64
    alpha_cap::Float64
    clip_to_coast::Symbol
    candidate_cache_limit::Int
    land_wkb::Vector{Cuchar}
    s2_stats
end

mutable struct S2ClipStats
    counts::Dict{Symbol,Atomic{Int}}
    timings_ns::Dict{Symbol,Atomic{Int64}}
end

const _S2_STAT_NAMES = (
    :interior_covers, :interior_miss, :regional_candidate,
    :regional_cache_hit, :regional_cache_miss, :regional_empty_candidates,
    :regional_cache_reset, :fallback_candidate_disabled, :fallback_global_extent,
    :fallback_buffer_failure, :fallback_query_failure, :fallback_index_failure,
    :fallback_too_many_candidates, :regional_empty_result,
    :fallback_local_clip_failure, :global_overlay,
)

const _S2_TIMING_NAMES = (
    :interior_total, :regional_index, :regional_buffer, :regional_tree_query,
    :regional_index_keys, :regional_cache_lookup, :regional_candidate_build,
    :regional_overlay, :regional_total, :global_total,
)

S2ClipStats() = S2ClipStats(
    Dict(name => Atomic{Int}(0) for name in _S2_STAT_NAMES),
    Dict(name => Atomic{Int64}(0) for name in _S2_TIMING_NAMES),
)

@inline function _s2_stat!(stats::S2ClipStats, name::Symbol)
    atomic_add!(stats.counts[name], 1)
    return nothing
end

@inline function _s2_time!(stats::S2ClipStats, name::Symbol, elapsed_ns::Integer)
    atomic_add!(stats.timings_ns[name], Int64(elapsed_ns))
    return nothing
end

s2_clip_stats(context::RangeContext) =
    Dict(name => counter[] for (name, counter) in context.s2_stats.counts)

s2_clip_timings(context::RangeContext) =
    Dict(name => counter[] / 1e9 for (name, counter) in context.s2_stats.timings_ns)

mutable struct WorkflowSpatialContext
    context
    forward
    inverse
    alpha_forward
    alpha_inverse
    geos::LG.GEOSContext
end

struct SpatialWorker
    transforms::WorkflowSpatialContext
    geos::LG.GEOSContext
    land_cea::LG.Geometry
    land_wgs84::LG.Geometry
    land_s2::Any
    land_prepared::Any
    land_cea_components::Vector{Any}
    land_wgs84_components::Vector{Any}
    land_cea_tree::LG.STRtree
    land_component_index::IdDict{Any,Int}
    land_candidate_cache::Dict{BitSet,Any}
    candidate_cache_limit::Int
    land_candidate_enabled::Bool
    s2_stats::S2ClipStats
    grid_geometries::Vector{Any}
    tree::LG.STRtree
    geometry_index::IdDict{Any,Int}
end

mutable struct DirectWorker
    transforms::WorkflowSpatialContext
    geos::LG.GEOSContext
    land_cea::LG.Geometry
    land_prepared::Any
    grid_geometries::Union{Nothing,Vector{Any}}
    grid_tree::Any
    candidate_geometries::Dict{Int,Any}
end

"""Normalized occurrence columns consumed by the bundled workflow."""
struct OccurrenceInput
    species::Vector{String}
    longitude::Vector{Float64}
    latitude::Vector{Float64}
    metadata::NamedTuple
    source_rows::Vector{Int}
end

Base.length(input::OccurrenceInput) = length(input.species)

struct SpeciesTask
    name::Symbol
    points::Matrix{Float64}
    metadata::NamedTuple
    source_rows::Vector{Int}
end

struct SpeciesRangeResult
    species::Symbol
    point_count::Int
    method::Symbol
    status::Symbol
    success::Bool
    metadata::NamedTuple
    source_rows::Vector{Int}
    grid_ids::Vector{Int}
    range_area_km2::Union{Missing,Float64}
    area_status::Symbol
    geometry
    error_message::Union{Nothing,String}
end

struct SpeciesRangesResult
    results::Vector{SpeciesRangeResult}
    elapsed_seconds::Float64
    s2_stats::Dict{Symbol,Int}
    s2_timings::Dict{Symbol,Float64}
end
