module Workflow

"""Optimized, bundled spatial workflow for multi-species ranges.

The workflow accepts named collections of `(longitude, latitude)` matrices and
uses only package-bundled spatial resources. It intentionally has no file
format or GBIF dependency.
"""

using Base.Threads

const _RB = parentmodule(@__MODULE__)
const GI = _RB.GeoInterface
const GO = _RB.GeometryOps
const LG = _RB.LibGEOS
const Proj = _RB.Proj
const S2 = _RB.S2Geography
const Tables = _RB.Tables

const EQUAL_AREA_CRS = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
const REGIONAL_LAND_MARGIN_M = 250_000.0
const REGIONAL_LAND_MAX_FRACTION = 0.75
const BUFFER_QUADRANT_SEGMENTS = 10
const _EARTH_WIDTH_M = 40_075_016.68557849
const _EARTH_HEIGHT_M = 20_003_931.45862545

export GridData, GridLocator, RangeContext, OccurrenceInput, SpeciesRangeResult,
       SpeciesRangesResult, clean_species_points, read_one_degree_map,
       bundled_land, new_range_context, prepare_workflow_context,
       standardize_occurrences, build_species_ranges,
       clip_hull_to_context, intersecting_grid_ids,
       s2_clip_stats, s2_clip_timings

include("types.jl")
include("resources.jl")
include("spatial.jl")
include("clipping.jl")
include("scheduler.jl")
include("results.jl")

end
