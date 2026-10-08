# Bundled range workflow

`RangeBuilder.Workflow` is the package-native multi-species range workflow.
Its primary input is a Tables.jl-compatible occurrence table. It also retains
a lower-level `NamedTuple`/dictionary matrix interface for compatibility. It
does not read files and does not depend on `ArchGDAL`, `CSV`, or `DataFrames`.

```julia
input = standardize_occurrences(
    df1;
    species=:Accepted_species,
    longitude=:mylongitude,
    latitude=:mylatitude,
    family=:Accepted_family,
    order=nothing,
)
```

`coordinate_digits` controls coordinate precision before de-duplication. It is
the number of decimal places retained for both longitude and latitude and
defaults to `2`; pass `1`, `3`, or another non-negative integer as needed. The
de-duplication key is `(species, rounded_longitude, rounded_latitude)`, and the
rounded coordinates are also the values used by the downstream range workflow.
Thus, duplicate rows after rounding are collapsed while `source_rows` retains
the original row index of the first retained record.

`input` is an `OccurrenceInput` with standardized `species`, `longitude`, and
`latitude` vectors, a `metadata` named tuple, and original `source_rows`. Invalid
coordinates and blank species names are removed during standardization, as are
duplicates after coordinate rounding; `source_rows` points back to the first
retained source row.

```julia
using RangeBuilder

context = prepare_workflow_context(
    buffer_m=10_000,
    clip_to_coast=:terrestrial,
    candidate_cache_limit=1024,
)

result = build_species_ranges(input, context; workers=8)
```

The context uses the bundled one-degree grid and Natural Earth 1:50m land.
Each worker owns its PROJ, GEOS, S2, and STRtree native state. Terrestrial
ranges first try the prepared CEA `covers` shortcut, then a conservative
regional land-component query and cached S2 geography, and finally the global
S2 overlay when the regional path is unsafe.

`candidate_cache_limit` bounds the number of cached regional candidate
geographies per worker. A value of `0` disables the cache. When the limit is
reached, the worker clears the cache before inserting the next candidate; this
controls native geometry memory without changing geometry semantics.

`SpeciesRangeResult.status` is currently `:success` or `:failed`. Species with
fewer than three valid unique points use the direct point/grid path. Every
species with at least three points, including collinear triples, is passed to
the R-compatible SHull/MCH recovery path.
Each result also carries the collapsed species-level `metadata` and exact
`source_rows` into the standardized input.

The compatibility entry point remains available when callers already have
species matrices:

```julia
points = Dict(:oak => [116.0 39.0; 116.5 39.0; 116.4 39.4])
result = build_species_ranges(points, context; workers=1)
```

On the Rosales reference input (`E:\Rosales\native_records.csv.gz`, 1,063,070
valid records and 9,263 species), the complete table-based run with 8 workers
finished with 9,263/9,263 successes in `build_ranges=101.818 s` and total
elapsed time `124.404 s`. Its point counts, grid counts, methods, areas, and
grid memberships matched `R/output/pluto_1deg_ranges/` exactly. The raw output
CSV SHA-256 is not used as the workflow contract because the new writer includes
`status` and omits the legacy `occurrence_sources` column; the geometry and
numeric result comparison is the relevant equivalence check.

For memory control, each worker owns its native PROJ/GEOS/S2/STRtree state and
the regional candidate cache is bounded by `candidate_cache_limit` (default
1024). Cache eviction only releases reusable native geometries; it does not
change the global-S2 fallback or result semantics.

The API intentionally leaves file I/O and tabular output to the application
layer. `SpeciesRangesResult.s2_stats` and `.s2_timings` provide diagnostics for
benchmarks; they are not a substitute for full-output equivalence validation.
