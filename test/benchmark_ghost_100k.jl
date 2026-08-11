using RangeBuilder
using DelaunayTriangulation
using Random

@eval DelaunayTriangulation begin
    function has_ghost_vertices(G::Graph{I}) where {I}
        return I(𝒢) in get_vertices(G)
    end
end

n = 100_000
points = rand(MersenneTwister(20260809), n, 2)
println("n = ", n, "; sentinel ghost check")
delvor(points)
GC.gc()
dv = @timed delvor(points)
println("delvor: ", round(dv.time; digits=4), " s; allocated: ",
        round(dv.bytes / 2.0^20; digits=2), " MiB; mesh rows: ", size(dv.value.mesh, 1))
ahull(dv.value; alpha=0.07)
GC.gc()
hull = @timed ahull(dv.value; alpha=0.07)
println("ahull cached: ", round(hull.time; digits=4), " s; arcs: ", size(hull.value.arcs, 1))
areaahull(hull.value)
GC.gc()
area = @timed areaahull(hull.value)
println("areaahull cached: ", round(area.time; digits=4), " s; area: ", area.value)
