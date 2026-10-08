function Base.show(io::IO, dv::DelVor)
    print(io, "DelVor(", size(dv.x, 1), " points, ", size(dv.mesh, 1), " mesh edges)")
end

function Base.show(io::IO, ::MIME"text/plain", dv::DelVor)
    println(io, "DelVor")
    println(io, "  x: ", size(dv.x), " Matrix{Float64}")
    println(io, "  mesh: ", size(dv.mesh), " Matrix{Float64}")
    print(io, "  triangulation: ", size(dv.x, 1), " points")
end

function Base.show(io::IO, sh::AShape)
    print(io, "AShape(", size(sh.edges, 1), " edges, alpha=", sh.alpha,
          ", length=", sh.length, ")")
end

function Base.show(io::IO, ::MIME"text/plain", sh::AShape)
    println(io, "AShape")
    println(io, "  edges: ", size(sh.edges), " Matrix{Float64}")
    println(io, "  length: ", sh.length)
    println(io, "  alpha: ", sh.alpha)
    println(io, "  alpha_extremes: ", sh.alpha_extremes)
    print(io, "  delvor: ", sh.delvor)
end

function Base.show(io::IO, h::AHull)
    print(io, "AHull(", size(h.arcs, 1), " arcs, alpha=", h.alpha,
          ", length=", h.length, ")")
end

function Base.show(io::IO, ::MIME"text/plain", h::AHull)
    println(io, "AHull")
    println(io, "  arcs: ", size(h.arcs), " Matrix{Float64}")
    println(io, "  xahull: ", size(h.xahull), " Matrix{Float64}")
    println(io, "  length: ", h.length)
    println(io, "  alpha: ", h.alpha)
    println(io, "  complement: ", size(h.complement), " Matrix{Float64}")
    print(io, "  ashape: ", h.ashape)
end

function Base.show(io::IO, ci::CircleIntersection)
    print(io, "CircleIntersection(n_cut=", ci.n_cut,
          ", theta1=", ci.theta1, ", theta2=", ci.theta2, ")")
end

function Base.show(io::IO, ::MIME"text/plain", ci::CircleIntersection)
    println(io, "CircleIntersection")
    println(io, "  n_cut: ", ci.n_cut)
    println(io, "  v1: ", ci.v1, ", theta1: ", ci.theta1)
    print(io, "  v2: ", ci.v2, ", theta2: ", ci.theta2)
end
