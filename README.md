# RangeBuilder.jl

RangeBuilder.jl is a high-performance Julia package for alpha-hull geometry, species range construction, and spatial richness analysis. Its alpha-hull kernel uses a pure-Julia SHull triangulation.

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/wyx619/RangeBuilder.jl")
```

## Quick start

```julia
using RangeBuilder

points = [
    116.0 39.0
    116.5 39.0
    116.4 39.4
    116.1 39.3
]

hull = ahull(points; alpha=0.7)
polygon = ah2polygon(hull)
inside = inahull(hull, [116.2 39.2; 118.0 40.0])

range = getDynamicAlphaHull(
    points;
    fraction=0.95,
    partCount=3,
    buff=1_000,
    clipToCoast=:terrestrial,
)

range.alpha
range.hull
```

Coordinates are supplied as an `n x 2` matrix in `(longitude, latitude)` order. `buff` is measured in metres. The result contains a GeoInterface-compatible polygon or multipolygon and the selected alpha label.

## Performance

RangeBuilder uses SHull triangulation for alpha-hull construction. Dynamic
alpha search reuses a successful non-jittered triangulation across candidates,
with recovery paths for duplicate points and degenerate geometry.

The range workflow is optimized around the remaining spatial hotspots:

- a CEA prepared-land `covers` check skips S2 overlay for fully terrestrial ranges;
- coastal ranges query a thread-local land-component `STRtree` before local S2 overlay, with automatic fallback to global S2 where required for safety;
- workers reuse PROJ, GEOS, S2, and grid spatial-index resources;
- unused global grid-centre construction has been removed.

On the Rosales reference input (1,063,070 records; 9,263 species), the bundled
Tables.jl workflow completed **9,263/9,263** species successfully with 8 Julia
workers: `build_ranges=101.818 s`, total elapsed time `124.404 s`.

## Dynamic range API

```julia
range = getDynamicAlphaHull(
    points;
    fraction=0.95,
    partCount=3,
    initialAlpha=2,
    alphaIncrement=1,
    alphaCap=400,
    buff=10_000,
    clipToCoast=:terrestrial,
)
```

`clipToCoast` accepts `:no`, `:terrestrial`, or `:aquatic`. The dynamic search proceeds within each species; independent species can be scheduled with Julia threads.

## S2 coast semantics

Accepted ranges are buffered in Equal Earth and then converted to WGS84 for the final geographic coast overlay. `S2Geography.jl` performs the terrestrial or aquatic intersection using Google S2 semantics, handling great-circle boundaries, antimeridian crossings and multi-component output.

The S2 overlay is deliberately limited to final coast clipping. The high-frequency alpha validity and occurrence-coverage checks use the lighter GeometryOps spherical predicates, while the dedicated spherical convex-hull implementation is used for the MCH fallback.

## Bundled spatial data

The package bundles Natural Earth 1:50m land data for offline coast clipping and land filtering. The one-degree grid used by the optional range-building workflow is also bundled as an internal JLD2 resource. Its path is resolved by the package; callers do not need to locate or pass the file.

## Multiple species and richness

The package also includes a bundled, file-free multi-species range workflow.
It accepts any Tables.jl-compatible occurrence table:

```julia
context = prepare_workflow_context()
input = standardize_occurrences(
    df1;
    species=:Accepted_species,
    longitude=:mylongitude,
    latitude=:mylatitude,
    family=:Accepted_family,
    order=nothing,
)
result = build_species_ranges(input, context; workers=8)
```

`standardize_occurrences` rounds longitude and latitude before de-duplication.
The default is two decimal places; set `coordinate_digits=1` (or any other
non-negative integer) to choose a different precision. De-duplication uses
`(species, rounded_longitude, rounded_latitude)`, and the rounded coordinates
are passed to the range calculation. `source_rows` records the original row
index of the first row retained for each rounded coordinate.

It uses the bundled one-degree grid and Natural Earth land data, with prepared
CEA land checks, regional S2 candidate caching, and safe global S2 fallback.
Each `SpeciesRangeResult` retains collapsed metadata and the source row indices
from the standardized input. The workflow does not perform file I/O; use the
application's table reader and writer around this API.
See [`src/workflow/README.md`](src/workflow/README.md) for the workflow API.

```julia
ranges = (
    species_a=ah2polygon(ahull(points; alpha=0.7)),
    species_b=ah2polygon(ahull(points .+ [0.5 0.5]; alpha=0.7)),
)

stack = rasterStackFromPolyList(ranges; resolution=0.2)
richness = speciesRichness(stack)
```

`rasterStackFromPolyList` returns a `Rasters.RasterStack`; `speciesRichness` reduces it to per-cell species counts. The package does not include a plotting backend.

## External plotting

Polygon outputs implement GeoInterface and can be plotted by the application that owns the display stack. For example, an application may install `CairoMakie` and `GeoMakie` and pass `range.hull` to `poly!`. Keeping plotting external avoids forcing a graphics backend on batch and server environments.
