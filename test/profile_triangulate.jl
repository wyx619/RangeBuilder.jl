using DelaunayTriangulation
using Random
using Profile

n = 20_000
points = rand(MersenneTwister(20260809 + n), n, 2)
pts = [(points[i, 1], points[i, 2]) for i in axes(points, 1)]
triangulate(pts; rng=MersenneTwister(7))
GC.gc()
Profile.clear()
@profile triangulate(pts; rng=MersenneTwister(7))
Profile.print(stdout; format=:flat, sortedby=:count, mincount=5)
