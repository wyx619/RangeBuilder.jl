# RangeBuilder.jl

A JuliaGeometry implementation of the computational API in R `alphahull`.
The geometry kernel is independent of R and uses `DelaunayTriangulation.jl`.
`RCall.jl` is only used by the reference regression test.

RangeBuilder.jl includes a runtime compatibility fast path for
`DelaunayTriangulation.jl 1.6.6`: ghost-vertex presence is checked by an
O(1) sentinel lookup instead of repeatedly scanning all graph vertices. No
upstream package modification is required.

## Performance

Performance is a primary design goal. RangeBuilder.jl is optimized for the
two expensive stages of alpha-hull range construction: Delaunay/Voronoi
construction and repeated candidate evaluation during dynamic alpha search.

- **One triangulation per dynamic search.** `getDynamicAlphaHull` builds the
  Delaunay/Voronoi structure once, then evaluates every candidate alpha against
  that cached structure. A search over *k* alpha values therefore performs one
  triangulation, rather than *k* triangulations.
- **A protected Delaunay fast path.** The `DelaunayTriangulation.jl 1.6.6`
  compatibility path replaces a repeated whole-graph ghost-vertex scan with an
  O(1) sentinel lookup. It prevents the severe scaling regression observed in
  large point sets without patching the upstream package.
- **Optimized alpha-hull assembly.** Arc-circle prefiltering, cached
  complement indices, and endpoint-indexed arc ordering remove avoidable work
  from the alpha-hull construction path. Coverage checks use prepared GEOS
  geometries; range buffering reuses coordinate transformations.

On the development machine, a warmed dynamic search over 10,000 points and
four candidate alphas completed in approximately **1.5–1.6 seconds** after
these changes, versus about **46.7 seconds** for the earlier implementation.
Timing is data- and hardware-dependent, but the structural improvement is
intentional: costly triangulation is no longer repeated for every alpha.

### Comparison with the original R packages

The following warmed median timings were measured on the same development
machine using deterministic inputs. They compare this package directly with
the installed original R `alphahull` and `rangeBuilder` packages—not with a
reimplementation. They are deliberately small enough to be rerun locally.

| Workload | Julia RangeBuilder.jl | Original R package | Relative result |
| --- | ---: | ---: | ---: |
| `ahull` + `areaahull`, 500 points | 0.0071 s | `alphahull`: 0.0600 s | 8.5× faster |
| `ahull` + `areaahull`, 2,000 points | 0.0318 s | `alphahull`: 0.2200 s | 6.9× faster |
| `ahull` + `areaahull`, 10,000 points | 0.2225 s | `alphahull`: 1.1300 s | 5.1× faster |
| Dynamic search, 500 points, four alpha candidates | 0.0453 s | `rangeBuilder`: 4.17 s | 92× faster |
| Dynamic search, 1,000 points, ten alpha candidates | 0.4474 s | `rangeBuilder`: 30.75 s | 68.7× faster |

The dynamic benchmarks use `fraction=0.95`, `partCount=3`, `buff=0`, and
disabled coast clipping so they isolate range construction. The 500-point row
searches `0.005:0.005:0.020`; the 1,000-point row searches
`0.005:0.005:0.050`. Both implementations chose the same fallback label,
`alphaMCH`, for each fixed point set. The large dynamic-search gap is
expected: the original R workflow repeatedly invokes `alphahull::ahull()` as
alpha changes, whereas RangeBuilder.jl reuses its single Delaunay/Voronoi
construction.

The benchmark scripts under `test/` measure both cached alpha-hull components
and end-to-end construction at fixed random seeds. They are designed for
revision-to-revision comparisons on the same machine.

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

## External plotting

RangeBuilder.jl intentionally has no plotting dependency and exports no
plotting API. This keeps the computational environment small and avoids
forcing a graphics backend on batch or server workflows. Its polygon outputs
implement the GeoInterface geometry protocol, and its richness output is a
`Rasters.Raster`; both can be plotted by packages installed in the calling
project.

For example, install a Makie backend and GeoMakie in an application that
needs maps:

```julia
import Pkg
Pkg.add(["CairoMakie", "GeoMakie"])
```

Then draw a computed range and its source occurrences without adding either
package to RangeBuilder.jl itself:

```julia
using CairoMakie, GeoMakie

range = getDynamicAlphaHull(points; buff=1_000, clipToCoast=:terrestrial)
fig = Figure()
ax = GeoAxis(fig[1, 1]; dest="+proj=eqearth")
poly!(ax, range.hull; color=(:steelblue, 0.35), strokecolor=:steelblue)
scatter!(ax, points[:, 1], points[:, 2]; color=:black, markersize=7)
fig
```

Use a different Makie backend (for example `GLMakie` or `WGLMakie`) when an
interactive desktop or browser display is needed. For a simple non-geographic
view of a richness raster, its values can be materialized for any plotting
library, for example `heatmap(Array(richness))` with Makie.

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

## Verification

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
julia --project=test test/rcall_reference.jl
julia --project=. test/benchmark_core.jl
julia --project=. test/benchmark_end_to_end.jl
Rscript test/benchmark_end_to_end.R
julia --project=test test/benchmark_dynamic_reference.jl
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

`benchmark_dynamic_reference.jl` compares `getDynamicAlphaHull` with the
installed original R `rangeBuilder` package through RCall. RCall is confined
to the test environment and remains outside the package's runtime
dependencies.
