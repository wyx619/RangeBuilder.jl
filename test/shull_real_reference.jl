# Optional real-data differential test. Run explicitly with:
# $env:RUN_REAL_SHULL_TEST='true'; julia --project=test test/shull_real_reference.jl
using RCall
using Test
using RangeBuilder

const REAL_OCCURRENCE_FILE = "E:/Rosales/native_records.csv.gz"

function r_real_species_points(species::String)
    isfile(REAL_OCCURRENCE_FILE) || error("occurrence file not found: $REAL_OCCURRENCE_FILE")
    occurrence_file = REAL_OCCURRENCE_FILE
    @rput occurrence_file species
    RCall.reval(
        "library(data.table); " *
        ".real_occurrences <- fread(occurrence_file, " *
        "select=c('Accepted_name', 'decimalLongitude', 'decimalLatitude'), " *
        "encoding='UTF-8', showProgress=FALSE); " *
        ".real_points <- unique(.real_occurrences[Accepted_name == species, " *
        ".(decimalLongitude, decimalLatitude)]); " *
        ".real_points <- .real_points[is.finite(decimalLongitude) & " *
        "is.finite(decimalLatitude) & decimalLongitude >= -180 & " *
        "decimalLongitude <= 180 & decimalLatitude >= -90 & decimalLatitude <= 90]",
    )
    return Matrix{Float64}(RCall.rcopy(RCall.reval("as.matrix(.real_points)")))
end

function r_real_delvor_mesh(points::Matrix{Float64})
    x = points[:, 1]
    y = points[:, 2]
    @rput x y
    RCall.reval("library(alphahull)")
    return Float64.(RCall.rcopy(RCall.reval("alphahull::delvor(cbind(x, y))\$mesh")))
end

function r_real_ahull_arcs(points::Matrix{Float64}, alpha::Float64)
    x = points[:, 1]
    y = points[:, 2]
    @rput x y alpha
    RCall.reval("library(alphahull)")
    return Float64.(RCall.rcopy(RCall.reval("alphahull::ahull(cbind(x, y), alpha=alpha)\$arcs")))
end

function r_real_dynamic_alpha(points::Matrix{Float64}; buff::Real=0.0)
    x = points[:, 1]
    y = points[:, 2]
    buff = Float64(buff)
    @rput x y buff
    RCall.reval("library(rangeBuilder)")
    return String(RCall.rcopy(RCall.reval(
        "rangeBuilder::getDynamicAlphaHull(cbind(x, y), buff=buff, clipToCoast='no')[[2]]",
    )))
end

if get(ENV, "RUN_REAL_SHULL_TEST", "false") == "true"
@testset "SHull real species references" begin
    for species in (
        "Alchemilla parcipila",
        "Ficus abscondita",
        "Rosa acicularis subsp. acicularis",
        "Potentilla elegans",
        "Prunus guanaiensis",
    )
        points = r_real_species_points(species)
        @test size(points, 1) >= 3
        julia_mesh = RangeBuilder.delvor(points; backend=:shull).mesh
        @test isapprox(julia_mesh, r_real_delvor_mesh(points); rtol=2e-6, atol=1e-10)
        julia_arcs = RangeBuilder.ahull(RangeBuilder.delvor(points; backend=:shull); alpha=3.0).arcs
        @test isapprox(julia_arcs, r_real_ahull_arcs(points, 3.0); rtol=2e-6, atol=1e-10)
    end

    for species in (
        "Alchemilla parcipila",
        "Ficus abscondita",
        "Rosa acicularis subsp. acicularis",
        "Potentilla biflora",
    )
        points = r_real_species_points(species)
        julia_dynamic = RangeBuilder.getDynamicAlphaHull(points; buff=0, clipToCoast=:no, backend=:shull)
        @test julia_dynamic.alpha == r_real_dynamic_alpha(points)
        julia_buffered_dynamic = RangeBuilder.getDynamicAlphaHull(
            points; buff=10_000, clipToCoast=:no, backend=:shull,
        )
        @test julia_buffered_dynamic.alpha == r_real_dynamic_alpha(points; buff=10_000)
    end
end
end
