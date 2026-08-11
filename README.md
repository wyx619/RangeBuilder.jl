# RangeBuilder.jl

A JuliaGeometry implementation of the computational API in R `alphahull`.
The geometry kernel is independent of R and uses `DelaunayTriangulation.jl`.
`RCall.jl` is only used by the reference regression test.

RangeBuilder.jl includes a runtime compatibility fast path for
`DelaunayTriangulation.jl 1.6.6`: ghost-vertex presence is checked by an
O(1) sentinel lookup instead of repeatedly scanning all graph vertices. No
upstream package modification is required.

## Use

```julia
using RangeBuilder
points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
hull = ahull(points; alpha=0.7)
inside = inahull(hull, [0.5 0.5; 3.0 3.0])

range = getDynamicAlphaHull(points; fraction=0.95, clipToCoast=:no)
filtered = filterByProximity(points, 20.0)

ranges = (species_a=ah2polygon(ahull(points; alpha=0.7)),
          species_b=ah2polygon(ahull(points .+ 0.5; alpha=0.7)))
range_stack = rasterStackFromPolyList(ranges; resolution=0.2)
richness = speciesRichness(range_stack)

batch = buildRanges(Dict("species_a" => points, "species_b" => points[1:2, :]);
                    clipToCoast=:no)
```

`getDynamicAlphaHull` returns `(hull, alpha)`, where `hull` is a
GeoInterface polygon/multipolygon and `alpha` records the selected alpha (or
`"alphaMCH"` for the convex-hull fallback). `buff` is measured in metres;
coast clipping uses the bundled Natural Earth 50m coastline data and works
offline.

`rasterStackFromPolyList` accepts a named tuple or dictionary of GeoInterface
polygons and returns a `Rasters.RasterStack`. Each layer contains `1` where a
cell center is inside the polygon and `missing` elsewhere. Use
`retainSmallRanges=false` to drop polygons with no covered cell centers, or
provide `extent=[xmin, xmax, ymin, ymax]` to override the automatic extent.
`speciesRichness` reduces that stack to one raster of per-cell species counts;
with its default `zeroToMissing=true`, cells with no species are `missing` as
in the R example workflow.

`filterByLand` returns `true` for occurrence coordinates on a Natural Earth
land polygon buffered by 2 km, `false` for ocean, and `missing` for missing
coordinates. The repository bundles Natural Earth v5.1.2's 50m land dataset,
so coastline clipping and land filtering work fully offline.

`buildRanges` is the conservative multi-species entry point. It returns
successfully modelled polygons in `ranges` and records species with fewer than
three unique points or degenerate geometry in `excluded`; these species are
not converted to inferred buffer ranges.

## Compatibility

The following computational functions preserve R's matrix schemas and
one-based point indices: `delvor`, `ashape`, `complement`, `ahull`, `dw`,
`inter`, `lengthahull`, `areaahulleval`, `areaahull`, `inahull`, `rotation`,
`anglesArc`, `koch`, and `rkoch`.

`arc` returns its sampled coordinates rather than drawing to a graphics
device. Similarly, `dw_track` and `ahull_track` return a vector of sampled
arc-coordinate matrices instead of R `ggplot2` layers. They accept an `rng`
keyword for reproducible sampling.

The Delaunay implementations use different internal traversal orders in R and
Julia. For non-degenerate inputs the resulting geometry is equivalent, though
the order of Delaunay rows and inserted alpha-hull endpoint ids can differ.

## Plotting

Plotting is provided by an optional Makie extension. Install and load a Makie
backend in the calling environment; it is not a runtime dependency of
RangeBuilder.jl. For static files, CairoMakie is a suitable backend:

```julia
using RangeBuilder
using CairoMakie

points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
hull = ahull(points; alpha=0.7)
fig = CairoMakie.plot(hull; do_shape=true, wlines=:vor, number=true)
save("ahull.png", fig)
```

`plot` and `plot!` are available through any loaded Makie backend for
`DelVor`, `AShape`, and `AHull`. The plotting keywords mirror R: `wlines` accepts `:none`, `:both`,
`:del`, or `:vor`; `wpoints`, `number`, `col`, `lwd`, `xlim`, and `ylim` are
also supported. `AHull` additionally accepts `do_shape`.

`plotdata` returns the backend-independent drawing commands used by the Makie
extension. It is useful when rendering with another Julia plotting package.

## Verification

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
julia --project=test test/rcall_reference.jl
julia --project=. test/benchmark_core.jl
julia --project=. test/benchmark_end_to_end.jl
Rscript test/benchmark_end_to_end.R
```

The RCall reference script compares against the installed `alphahull` R
package, not the repository source files. RCall is isolated in
`test/Project.toml`, so it is not a runtime dependency of RangeBuilder.jl.

`benchmark_core.jl` warms up Julia and reports median time and allocation for
`delvor`, `ashape`, `complement`, `ahull`, and `dw` at fixed random seeds.
The benchmark is intended for comparing revisions on the same machine, not
for cross-machine timing claims.

`benchmark_end_to_end.jl` and `benchmark_end_to_end.R` use the same
deterministic point sets for a warm end-to-end comparison of `ahull` plus
area evaluation.
