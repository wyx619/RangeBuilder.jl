using DataFrames
using Random
using Statistics

include(joinpath(@__DIR__, "..", "R", "build_species_alpha_hull_1deg_ranges_gbif.jl"))
using .RangeBuilderWorkflow

const GRID_PATH = raw"C:\Users\MSI\GitHub\alphahull\R\1d\1.shp"
const OCCURRENCE_PATH = raw"E:\Rosales\native_records.csv.gz"
const SAMPLE_SIZE = 93
const SAMPLE_SEED = 20260813

grid = read_range_grid(GRID_PATH)
context = new_range_context(
    grid;
    buffer_m=10_000,
    fraction=0.95,
    part_count=3,
    initial_alpha=2,
    alpha_increment=1,
    alpha_cap=400,
    clip_to_coast=:terrestrial,
)
points = read_occurrence_file(OCCURRENCE_PATH)
species = sort!(unique(points.species))
selection = sort!(randperm(MersenneTwister(SAMPLE_SEED), length(species))[1:SAMPLE_SIZE])
sample = filter(:species => in(Set(species[selection])), points)
counts = combine(groupby(sample, :species), nrow => :n)

# Compile every code path and construct per-task native spatial contexts.
build_family_ranges(sample, context; workers=8)
GC.gc()
single_time = @elapsed single = build_family_ranges(sample, context; workers=1)
GC.gc()
parallel_time = @elapsed parallel = build_family_ranges(sample, context; workers=8)

println("sample records: $(nrow(sample)); species: $(nrow(counts))")
println("point count (min / median / max): $(minimum(counts.n)) / $(median(counts.n)) / $(maximum(counts.n))")
println("workers=1: $(round(single_time; digits=3)) s")
println("workers=8: $(round(parallel_time; digits=3)) s")
println("speedup: $(round(single_time / parallel_time; digits=3))x")

single_keys = Set(zip(single.grids.species, single.grids.grid_id))
parallel_keys = Set(zip(parallel.grids.species, parallel.grids.grid_id))
single_area = Dict(row.species => row.range_area_km2 for row in eachrow(single.summary))
parallel_area = Dict(row.species => row.range_area_km2 for row in eachrow(parallel.summary))
single_success = Dict(row.species => row.success for row in eachrow(single.summary))
parallel_success = Dict(row.species => row.success for row in eachrow(parallel.summary))
status_mismatches = count(species -> single_success[species] != parallel_success[species], keys(single_success))
area_mismatches = count(keys(single_area)) do species
    left, right = single_area[species], parallel_area[species]
    ismissing(left) || ismissing(right) ? !ismissing(left) || !ismissing(right) :
        !isapprox(left, right; rtol=1e-10, atol=1e-6)
end
println("grid keys (single / parallel): $(length(single_keys)) / $(length(parallel_keys))")
println("grid keys only in single / parallel: $(length(setdiff(single_keys, parallel_keys))) / $(length(setdiff(parallel_keys, single_keys)))")
println("success-status mismatches: $status_mismatches")
println("area mismatches: $area_mismatches")
