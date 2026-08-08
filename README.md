# AlphaHull.jl

A JuliaGeometry implementation of the computational API in R `alphahull`.
The geometry kernel is independent of R and uses `DelaunayTriangulation.jl`.
`RCall.jl` is only used by the reference regression test.

## Use

```julia
using AlphaHull
points = [0.0 0.0; 2.0 0.0; 2.0 1.0; 0.0 2.0; 0.6 0.7]
hull = ahull(points; alpha=0.7)
inside = inahull(hull, [0.5 0.5; 3.0 3.0])
```

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
```

The RCall reference script compares against the installed `alphahull` R
package, not the repository source files. RCall is isolated in
`test/Project.toml`, so it is not a runtime dependency of AlphaHull.jl.

`benchmark_core.jl` warms up Julia and reports median time and allocation for
`delvor`, `ashape`, `complement`, `ahull`, and `dw` at fixed random seeds.
The benchmark is intended for comparing revisions on the same machine, not
for cross-machine timing claims.

`benchmark_end_to_end.jl` and `benchmark_end_to_end.R` use the same
deterministic point sets for a warm end-to-end comparison of `ahull` plus
area evaluation.
