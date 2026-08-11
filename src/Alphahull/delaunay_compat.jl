# DelaunayTriangulation 1.6.6 scans every graph vertex in this predicate.
# GhostVertex is always present when ghost vertices exist, so a hash lookup is
# equivalent and avoids an O(n) scan during every point-location step.
function _install_delaunay_ghost_fastpath!()
    method = which(DelaunayTriangulation.has_ghost_vertices,
                   Tuple{DelaunayTriangulation.Graph{Int}})
    Base.delete_method(method)
    @eval DelaunayTriangulation begin
        function has_ghost_vertices(G::Graph{I}) where {I}
            return I(GhostVertex) in get_vertices(G)
        end
    end
end

function __init__()
    _install_delaunay_ghost_fastpath!()
end
