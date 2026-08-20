# RangeBuilder.jl

RangeBuilder.jl is a high-performance JuliaGeometry implementation of the computational APIs from R `alphahull` and `rangeBuilder`. The geometry kernel is independent of R and does not require `RCall` at runtime.

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/wyx619/RangeBuilder.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path="C:/path/to/RangeBuilder.jl")
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

The main performance work is concentrated in Delaunay/Voronoi construction and dynamic alpha evaluation:

- one Delaunay structure is reused across candidate alphas in the cached backend;
- the `DelaunayTriangulation.jl 1.6.6` compatibility path uses an O(1) ghost-vertex sentinel lookup;
- alpha-hull assembly uses arc prefilters, cached complement indices and endpoint-indexed ordering;
- arc intersection candidates use sparse bounding-box adjacency rather than a dense `a x a` matrix;
- prepared GEOS predicates accelerate repeated coverage and grid-boundary checks without changing the boundary predicate.

The `:shull` backend rebuilds each candidate when exact R `rangeBuilder` call order is required. It is slower than the cached backend but preserves candidate-level retry and RNG semantics.

The development benchmark on 10,000 points and four candidate alphas is approximately 1.5-1.6 seconds after warm-up, compared with approximately 46.7 seconds for the earlier implementation. Timings depend on hardware and input geometry.

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
    backend=:shull,
)
```

`clipToCoast` accepts `:no`, `:terrestrial`, or `:aquatic`. The dynamic search itself remains ordered within one species; independent species can be scheduled by an application with Julia threads. For strict R-compatible retries, pass a matching `RSeed` RNG.

## S2 coast semantics

Accepted ranges are buffered in Equal Earth and then converted to WGS84 for the final geographic coast overlay. `S2Geography.jl` performs the terrestrial or aquatic intersection using Google S2 semantics. This handles great-circle boundaries, antimeridian crossings and multi-component output without R or `RCall`.

The S2 overlay is deliberately limited to final coast clipping. The high-frequency alpha validity and occurrence-coverage checks use the lighter GeometryOps spherical predicates, while the dedicated spherical convex-hull implementation is used for the MCH fallback.

## Bundled spatial data

The package bundles Natural Earth 1:50m land data for offline coast clipping and land filtering. The one-degree grid used by the optional range-building workflow is also bundled as an internal JLD2 resource. Its path is resolved by the package; callers do not need to locate or pass the file.

## Multiple species and richness

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

## Compatibility

The computational functions preserve R-style matrix schemas and one-based point indices for `delvor`, `ashape`, `complement`, `ahull`, `dw`, `inter`, `lengthahull`, `areaahulleval`, `areaahull`, `inahull`, `rotation`, `anglesArc`, `koch` and `rkoch`.

For ordinary inputs, R and Julia produce equivalent geometry. Cocircular or near-degenerate point sets can have multiple valid Delaunay triangulations, so internal triangle, arc and row order is not part of the compatibility contract. The accepted geometry and range semantics remain the target.

## Development checks

```julia
using Pkg
Pkg.test()
```


