using DataFrames
using Random
using Statistics

include(joinpath(@__DIR__, "..", "R", "build_species_alpha_hull_1deg_ranges_gbif.jl"))
using .RangeBuilderWorkflow

grid = read_range_grid(raw"C:\Users\MSI\GitHub\alphahull\R\1d\1.shp")
context = new_range_context(grid; buffer_m=10_000, fraction=0.95, part_count=3, initial_alpha=2, alpha_increment=1, alpha_cap=400, clip_to_coast=:terrestrial)
points = read_occurrence_file(raw"E:\Rosales\native_records.csv.gz")
species = sort!(unique(points.species))
selection = sort!(randperm(MersenneTwister(20260813), length(species))[1:93])
sample = filter(:species => in(Set(species[selection])), points)
observations = unique(select(copy(sample), RangeBuilderWorkflow.REQUIRED_COLUMNS), [:species, :longitude, :latitude])
tasks = sort!([DataFrame(group) for group in groupby(observations, :species) if nrow(group) >= context.alpha_min_points]; by=task -> task.species[1])

RangeBuilderWorkflow.map_alpha_hulls(tasks, context; workers=8)
GC.gc()
t1 = @elapsed alpha1 = RangeBuilderWorkflow.map_alpha_hulls(tasks, context; workers=1)
GC.gc()
t8 = @elapsed alpha8 = RangeBuilderWorkflow.map_alpha_hulls(tasks, context; workers=8)

GC.gc()
p1 = @elapsed single = RangeBuilderWorkflow.map_alpha_results(tasks, alpha1, context; workers=1)
GC.gc()
p8 = @elapsed parallel = RangeBuilderWorkflow.map_alpha_results(tasks, alpha8, context; workers=8)

single_keys = Set((row.species, row.grid_id) for row in eachrow(reduce(vcat, getproperty.(single, :grids))))
parallel_keys = Set((row.species, row.grid_id) for row in eachrow(reduce(vcat, getproperty.(parallel, :grids))))
single_area = Dict(result.summary.species[1] => result.summary.range_area_km2[1] for result in single)
parallel_area = Dict(result.summary.species[1] => result.summary.range_area_km2[1] for result in parallel)
area_mismatches = count(keys(single_area)) do species
    left, right = single_area[species], parallel_area[species]
    ismissing(left) || ismissing(right) ? !ismissing(left) || !ismissing(right) :
        !isapprox(left, right; rtol=1e-10, atol=1e-6)
end

println("alpha tasks: $(length(tasks)); records: $(nrow(sample)); max points/species: $(maximum(nrow.(tasks)))")
println("alpha search workers=1: $(round(t1; digits=3)) s")
println("alpha search workers=8: $(round(t8; digits=3)) s; speedup $(round(t1/t8; digits=3))x")
println("clip/grid postprocess workers=1: $(round(p1; digits=3)) s")
println("parallel clip/grid postprocess: $(round(p8; digits=3)) s; speedup $(round(p1/p8; digits=3))x")
println("grid keys only in single / parallel: $(length(setdiff(single_keys, parallel_keys))) / $(length(setdiff(parallel_keys, single_keys)))")
println("area mismatches: $area_mismatches")
