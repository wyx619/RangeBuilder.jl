"""Build the Float32 point representation used by `interp`'s SHull backend.

`interp::tri.mesh()` stores its working coordinates in C++ `float` fields
before sorting, seed selection, visibility tests, and edge flips. This helper
captures that topology precision separately from the Float64 coordinates used
by the alpha-hull output calculations.
"""
function _shull_topology_points(xy::AbstractMatrix)
    return [(Float32(xy[index, 1]), Float32(xy[index, 2]))
            for index in axes(xy, 1)]
end

"""Return the first pair that collapses to one SHull Float32 coordinate."""
function _shull_float32_duplicate_pair(xy::AbstractMatrix)
    seen = Dict{Tuple{Float32, Float32}, Int}()
    sizehint!(seen, size(xy, 1))
    for index in axes(xy, 1)
        key = (Float32(xy[index, 1]), Float32(xy[index, 2]))
        first_index = get(seen, key, 0)
        first_index == 0 || return (first_index, index)
        seen[key] = index
    end
    return nothing
end

"""Experimental Float32 topology proxy for the first SHull porting stage.

This intentionally does not replace `delvor`. It lets the test suite compare
the effect of SHull's Float32 topology precision against the locally installed
R `interp::tri.mesh()` before porting SHull's ordering and flip state machine.
"""
function _shull_topology_prototype(xy::AbstractMatrix)
    return triangulate(_shull_topology_points(xy); randomise=false)
end

mutable struct _SHullPoint
    id::Int
    trid::Int
    r::Float32
    c::Float32
    tr::Float32
    tc::Float32
    ro::Float32
end

_SHullPoint(id::Int, r::Float32, c::Float32) =
    _SHullPoint(id, 0, r, c, 0f0, 0f0, 0f0)

@inline function _copy_shull_point(point::_SHullPoint)
    return _SHullPoint(point.id, point.trid, point.r, point.c, point.tr, point.tc, point.ro)
end

mutable struct _SHullTriad
    a::Int
    b::Int
    c::Int
    ab::Int
    bc::Int
    ac::Int
    ro::Float32
    r::Float32
    ccentre::Float32
end

_SHullTriad(a::Int, b::Int, c::Int) =
    _SHullTriad(a, b, c, 0, 0, 0, -1f0, 0f0, 0f0)

"""Return SHull's Float32 circumcentre and squared radius.

This is a direct translation of `interp/src/s_hull_pro.cpp:circle_cent2`.
The Float64 coordinates retained by `delvor` are deliberately not used here:
SHull uses its `float` working representation for seed selection and topology.
"""
@inline function _shull_circle_cent2(r1::Float32, c1::Float32,
                                     r2::Float32, c2::Float32,
                                     r3::Float32, c3::Float32)
    a1 = (r1 + r2) / 2f0
    a2 = (c1 + c2) / 2f0
    b1 = (r3 + r2) / 2f0
    b2 = (c3 + c2) / 2f0
    e2 = r1 - r2
    e1 = -c1 + c2
    q2 = r3 - r2
    q1 = -c3 + c2

    e1 * -q2 + e2 * q1 == 0f0 && return (0f0, 0f0, -1f0)
    beta = (-e2 * (b1 - a1) + e1 * (b2 - a2)) / (e2 * q1 - e1 * q2)
    r = b1 + q1 * beta
    c = b2 + q2 * beta
    return (r, c, (r1 - r)^2 + (c1 - c)^2)
end

@inline _shull_less(left::_SHullPoint, right::_SHullPoint) =
    left.ro == right.ro ? (left.r == right.r ? left.c < right.c : left.r < right.r) :
    left.ro < right.ro

function _shull_r_jitter_amount(values::AbstractVector{<:Real}; factor::Real=0.001)
    finite_values = Float64[value for value in values if isfinite(value)]
    isempty(finite_values) && return 0.0
    lower, upper = extrema(finite_values)
    span = upper - lower
    span == 0 && (span = abs(lower))
    span == 0 && (span = 1.0)
    digits = 3 - floor(Int, log10(span))
    rounded = sort!(unique!(round.(finite_values; digits)))
    spacing = length(rounded) > 1 ? minimum(diff(rounded)) :
              (rounded[1] != 0 ? rounded[1] / 10 : span / 10)
    return Float64(factor) * abs(spacing) / 5
end

@inline _shull_uniform(rng::Random.AbstractRNG) = rand(rng)
@inline _shull_uniform(rng::RSeed.AbstractRUniformRNG) = RSeed.r_unif(rng)

function _shull_jittered_points(xy::AbstractMatrix;
                                rng=Random.MersenneTwister(0))
    # `interp::tri.mesh()` retries SHull failures with `jitter(x, 0.001)`.
    # As in R, x is perturbed first and y second. The deterministic local RNG
    # keeps the Julia backend reproducible while preserving that amplitude.
    jittered = Matrix{Float64}(xy)
    for column in axes(jittered, 2)
        amount = _shull_r_jitter_amount(view(jittered, :, column))
        @inbounds for row in axes(jittered, 1)
            jittered[row, column] += amount * (2 * _shull_uniform(rng) - 1)
        end
    end
    return jittered
end

@inline function _shull_requires_jitter(error)
    error isa ErrorException || return false
    message = sprint(showerror, error)
    return occursin("no visible hull facet", message) ||
           occursin("flip queue grew unexpectedly", message)
end

"""Build the pre-flip SHull triads used by `interp::tri.mesh()`.

This is the second, internal porting stage. It reproduces SHull's Float32
ordering, seed selection, visible-hull insertion, and triad creation exactly
enough to expose the deterministic state before Cline-Renka flips. It is not
used by `delvor`; the following flip stage must first prove ordered agreement
with R's `tri.mesh()`.
"""
function _shull_incremental_triads(xy::AbstractMatrix)
    size(xy, 2) == 2 || throw(ArgumentError("xy must have two columns"))
    nump = size(xy, 1)
    nump >= 3 || throw(ArgumentError("at least three points are required"))

    points = [_SHullPoint(index, Float32(xy[index, 1]), Float32(xy[index, 2]))
              for index in axes(xy, 1)]
    origin_r, origin_c = points[1].r, points[1].c
    for point in points
        dr = point.r - origin_r
        dc = point.c - origin_c
        point.ro = dr * dr + dc * dc
    end
    sort!(points; lt=_shull_less)

    r1, c1 = points[1].r, points[1].c
    r2, c2 = points[2].r, points[2].c
    middle = 0
    minimum_radius = 9f20
    centre_r = 0f0
    centre_c = 0f0
    index = 3
    while index <= nump
        r, c, radius = _shull_circle_cent2(r1, c1, r2, c2, points[index].r, points[index].c)
        if radius < minimum_radius && radius > 0f0
            middle = index
            minimum_radius = radius
            centre_r = r
            centre_c = c
        elseif minimum_radius * 4f0 < points[index].ro
            break
        end
        index += 1
    end
    middle > 0 || throw(ArgumentError("SHull could not select a non-collinear seed triangle"))

    point0 = _copy_shull_point(points[1])
    point1 = _copy_shull_point(points[2])
    point2 = _copy_shull_point(points[middle])
    deleteat!(points, middle)
    deleteat!(points, 1)
    deleteat!(points, 1)
    for point in points
        dr = point.r - centre_r
        dc = point.c - centre_c
        point.ro = dr * dr + dc * dc
    end
    sort!(points; lt=_shull_less)
    prepend!(points, [_copy_shull_point(point0), _copy_shull_point(point1), _copy_shull_point(point2)])

    centroid_r = (points[1].r + points[2].r + points[3].r) / 3f0
    centroid_c = (points[1].c + points[2].c + points[3].c) / 3f0
    dr0 = points[1].r - centroid_r
    dc0 = points[1].c - centroid_c
    tr01 = points[2].r - points[1].r
    tc01 = points[2].c - points[1].c

    hull = _SHullPoint[]
    triads = _SHullTriad[]
    if -tr01 * dc0 + tc01 * dr0 < 0f0
        point0.tr, point0.tc, point0.trid = point1.r - point0.r, point1.c - point0.c, 1
        point1.tr, point1.tc, point1.trid = point2.r - point1.r, point2.c - point1.c, 1
        point2.tr, point2.tc, point2.trid = point0.r - point2.r, point0.c - point2.c, 1
        append!(hull, [_copy_shull_point(point0), _copy_shull_point(point1), _copy_shull_point(point2)])
        triad = _SHullTriad(point0.id, point1.id, point2.id)
        triad.ro, triad.r, triad.ccentre = minimum_radius, centre_r, centre_c
        push!(triads, triad)
    else
        point0.tr, point0.tc, point0.trid = point2.r - point0.r, point2.c - point0.c, 1
        point2.tr, point2.tc, point2.trid = point1.r - point2.r, point1.c - point2.c, 1
        point1.tr, point1.tc, point1.trid = point0.r - point1.r, point0.c - point1.c, 1
        append!(hull, [_copy_shull_point(point0), _copy_shull_point(point2), _copy_shull_point(point1)])
        triad = _SHullTriad(point0.id, point2.id, point1.id)
        triad.ro, triad.r, triad.ccentre = minimum_radius, centre_r, centre_c
        push!(triads, triad)
    end

    for point_index in 4:nump
        candidate = _copy_shull_point(points[point_index])
        rx, cx = candidate.r, candidate.c
        numh = length(hull)
        numh > 0 || throw(ErrorException("SHull hull unexpectedly became empty"))
        dr = rx - hull[1].r
        dc = cx - hull[1].c
        visible = -dc * hull[1].tr + dr * hull[1].tc
        point_ids = Int[]
        triad_ids = Int[]
        hull_index = 0

        if visible < 0f0
            hull_index = 1
            dr = rx - hull[end].r
            dc = cx - hull[end].c
            visible = -dc * hull[end].tr + dr * hull[end].tc
            if visible < 0f0
                push!(point_ids, hull[end].id)
                push!(triad_ids, hull[end].trid)
                h = 1
                while h <= numh - 1
                    dr = rx - hull[h].r
                    dc = cx - hull[h].c
                    visible = -dc * hull[h].tr + dr * hull[h].tc
                    push!(point_ids, hull[h].id)
                    push!(triad_ids, hull[h].trid)
                    if visible < 0f0
                        deleteat!(hull, h)
                        numh -= 1
                    else
                        candidate.tr = hull[h].r - candidate.r
                        candidate.tc = hull[h].c - candidate.c
                        insert!(hull, 1, candidate)
                        numh += 1
                        break
                    end
                end
                h = numh - 1
                while h > 1
                    dr = rx - hull[h].r
                    dc = cx - hull[h].c
                    visible = -dc * hull[h].tr + dr * hull[h].tc
                    if visible < 0f0
                        insert!(point_ids, 1, hull[h].id)
                        insert!(triad_ids, 1, hull[h].trid)
                        deleteat!(hull, h + 1)
                    else
                        hull[end].tr = -hull[end].r + candidate.r
                        hull[end].tc = -hull[end].c + candidate.c
                        break
                    end
                    h -= 1
                end
            else
                hull_index = 2
                push!(triad_ids, hull[1].trid)
                push!(point_ids, hull[1].id)
                h = 2
                while h <= numh
                    dr = rx - hull[h].r
                    dc = cx - hull[h].c
                    visible = -dc * hull[h].tr + dr * hull[h].tc
                    push!(point_ids, hull[h].id)
                    push!(triad_ids, hull[h].trid)
                    if visible < 0f0
                        deleteat!(hull, h)
                        numh -= 1
                    else
                        candidate.tr = hull[h].r - candidate.r
                        candidate.tc = hull[h].c - candidate.c
                        hull[h - 1].tr = candidate.r - hull[h - 1].r
                        hull[h - 1].tc = candidate.c - hull[h - 1].c
                        insert!(hull, h, candidate)
                        break
                    end
                end
            end
        else
            # These indices retain the C++ zero-based convention because the
            # wrap-around branch uses `numh` as an exclusive end sentinel.
            first_visible = -1
            first_invisible = numh
            for h in 1:numh - 1
                dr = rx - hull[h + 1].r
                dc = cx - hull[h + 1].c
                visible = -dc * hull[h + 1].tr + dr * hull[h + 1].tc
                if visible < 0f0
                    first_visible < 0 && (first_visible = h)
                elseif first_visible > 0
                    first_invisible = h
                    break
                end
            end
            first_visible >= 0 || throw(ErrorException(
                "SHull found no visible hull facet at insertion $point_index",
            ))
            if first_invisible < numh
                for edge in first_visible:first_invisible
                    push!(point_ids, hull[edge + 1].id)
                    push!(triad_ids, hull[edge + 1].trid)
                end
            else
                for edge in first_visible:first_invisible - 1
                    push!(point_ids, hull[edge + 1].id)
                    push!(triad_ids, hull[edge + 1].trid)
                end
                push!(point_ids, hull[1].id)
            end
            if first_visible < first_invisible - 1
                deleteat!(hull, first_visible + 2:first_invisible)
            end
            if first_invisible == numh
                candidate.tr = hull[1].r - candidate.r
                candidate.tc = hull[1].c - candidate.c
            else
                candidate.tr = hull[first_visible + 2].r - candidate.r
                candidate.tc = hull[first_visible + 2].c - candidate.c
            end
            hull[first_visible + 1].tr = candidate.r - hull[first_visible + 1].r
            hull[first_visible + 1].tc = candidate.c - hull[first_visible + 1].c
            insert!(hull, first_visible + 2, candidate)
            hull_index = first_visible + 2
        end

        triangle_count = length(triads)
        edge_count = length(point_ids) - 1
        edge_count >= 1 || throw(ErrorException(
            "SHull produced an empty visible chain at insertion $point_index " *
            "(hull size $(length(hull)), initial visibility $visible)",
        ))
        if edge_count == 1
            triad = _SHullTriad(candidate.id, point_ids[1], point_ids[2])
            triad.bc = triad_ids[1]
            triad.ab = 0
            triad.ac = 0
            previous = triads[triad_ids[1]]
            if (triad.b == previous.a && triad.c == previous.b) ||
               (triad.b == previous.b && triad.c == previous.a)
                previous.ab = triangle_count + 1
            elseif (triad.b == previous.a && triad.c == previous.c) ||
                   (triad.b == previous.c && triad.c == previous.a)
                previous.ac = triangle_count + 1
            elseif (triad.b == previous.b && triad.c == previous.c) ||
                   (triad.b == previous.c && triad.c == previous.b)
                previous.bc = triangle_count + 1
            end
            hull[hull_index].trid = triangle_count + 1
            if hull_index > 1
                hull[hull_index - 1].trid = triangle_count + 1
            else
                hull[end].trid = triangle_count + 1
            end
            push!(triads, triad)
        else
            for edge in 1:edge_count
                triad = _SHullTriad(candidate.id, point_ids[edge], point_ids[edge + 1])
                triad.bc = triad_ids[edge]
                triad.ab = edge == 1 ? 0 : triangle_count
                triad.ac = triangle_count + 2
                previous = triads[triad_ids[edge]]
                if (triad.b == previous.a && triad.c == previous.b) ||
                   (triad.b == previous.b && triad.c == previous.a)
                    previous.ab = triangle_count + 1
                elseif (triad.b == previous.a && triad.c == previous.c) ||
                       (triad.b == previous.c && triad.c == previous.a)
                    previous.ac = triangle_count + 1
                elseif (triad.b == previous.b && triad.c == previous.c) ||
                       (triad.b == previous.c && triad.c == previous.b)
                    previous.bc = triangle_count + 1
                end
                push!(triads, triad)
                triangle_count += 1
            end
            triads[end].ac = 0
            hull[hull_index].trid = triangle_count
            if hull_index > 1
                hull[hull_index - 1].trid = length(triads) - edge_count + 1
            else
                hull[end].trid = length(triads) - edge_count + 1
            end
        end
    end
    return (; points, triads, hull)
end

@inline function _shull_triad(a::Int, b::Int, c::Int, ab::Int, bc::Int, ac::Int)
    return _SHullTriad(a, b, c, ab, bc, ac, -1f0, 0f0, 0f0)
end

@inline function _shull_cline_renka(a::_SHullPoint, b::_SHullPoint,
                                    c::_SHullPoint, d::_SHullPoint)
    v1r, v1c = b.r - a.r, b.c - a.c
    v2r, v2c = c.r - a.r, c.c - a.c
    v3r, v3c = b.r - d.r, b.c - d.c
    v4r, v4c = c.r - d.r, c.c - d.c
    cosine_a = v1r * v2r + v1c * v2c
    cosine_d = v3r * v4r + v3c * v4c
    cosine_a < 0f0 && cosine_d < 0f0 && return -1
    cosine_a > 0f0 && cosine_d > 0f0 && return 1
    sine_a = abs(v1r * v2c - v1c * v2r)
    sine_d = abs(v3r * v4c - v3c * v4r)
    return cosine_a * sine_d + sine_a * cosine_d < 0f0 ? -1 : 1
end

function _shull_replace_neighbour!(triad::_SHullTriad, old::Int, new::Int)
    triad.ab == old && (triad.ab = new; return)
    triad.bc == old && (triad.bc = new; return)
    triad.ac == old && (triad.ac = new; return)
    return
end

function _shull_bc_links(triad::_SHullTriad, neighbour::_SHullTriad, index::Int)
    if neighbour.ab == index
        links = triad.b == neighbour.a ?
                (neighbour.ac, neighbour.bc) : (neighbour.bc, neighbour.ac)
        return neighbour.c, links[1], links[2]
    elseif neighbour.ac == index
        links = triad.b == neighbour.a ?
                (neighbour.ab, neighbour.bc) : (neighbour.bc, neighbour.ab)
        return neighbour.b, links[1], links[2]
    elseif neighbour.bc == index
        links = triad.b == neighbour.b ?
                (neighbour.ab, neighbour.ac) : (neighbour.ac, neighbour.ab)
        return neighbour.a, links[1], links[2]
    end
    throw(ErrorException("SHull triad adjacency is inconsistent"))
end

function _shull_a_links(triad::_SHullTriad, neighbour::_SHullTriad, index::Int)
    if neighbour.ab == index
        links = triad.a == neighbour.a ?
                (neighbour.ac, neighbour.bc) : (neighbour.bc, neighbour.ac)
        return neighbour.c, links[1], links[2]
    elseif neighbour.ac == index
        links = triad.a == neighbour.a ?
                (neighbour.ab, neighbour.bc) : (neighbour.bc, neighbour.ab)
        return neighbour.b, links[1], links[2]
    elseif neighbour.bc == index
        links = triad.a == neighbour.b ?
                (neighbour.ab, neighbour.ac) : (neighbour.ac, neighbour.ab)
        return neighbour.a, links[1], links[2]
    end
    throw(ErrorException("SHull triad adjacency is inconsistent"))
end

"""Attempt one source-order SHull Cline-Renka edge flip.

`edge` follows the `Triad` fields from `s_hull_pro.cpp`: 1 is `bc`, 2 is
`ab`, and 3 is `ac`. Triad ids are one-based here; zero remains the boundary
sentinel used by the original C++ implementation.
"""
function _shull_try_flip!(points::Vector{_SHullPoint}, triads::Vector{_SHullTriad},
                          index::Int, edge::Int; require_distinct::Bool=true)
    triad = triads[index]
    neighbour_index = edge == 1 ? triad.bc : edge == 2 ? triad.ab : triad.ac
    neighbour_index == 0 && return 0
    neighbour = triads[neighbour_index]
    opposite, link3, link4 = edge == 1 ? _shull_bc_links(triad, neighbour, index) :
                             _shull_a_links(triad, neighbour, index)
    if edge == 1
        link1, link2 = triad.ab, triad.ac
        predicate = _shull_cline_renka(
            points[triad.a], points[triad.b], points[triad.c], points[opposite],
        )
        first = _shull_triad(triad.a, triad.b, opposite, link1, link3, neighbour_index)
        second = _shull_triad(triad.a, triad.c, opposite, link2, link4, index)
    elseif edge == 2
        link1, link2 = triad.ac, triad.bc
        predicate = _shull_cline_renka(
            points[triad.c], points[triad.b], points[triad.a], points[opposite],
        )
        first = _shull_triad(triad.c, triad.a, opposite, link1, link3, neighbour_index)
        second = _shull_triad(triad.c, triad.b, opposite, link2, link4, index)
    else
        link1, link2 = triad.ab, triad.bc
        predicate = _shull_cline_renka(
            points[triad.b], points[triad.a], points[triad.c], points[opposite],
        )
        first = _shull_triad(triad.b, triad.a, opposite, link1, link3, neighbour_index)
        second = _shull_triad(triad.b, triad.c, opposite, link2, link4, index)
    end
    predicate < 0 || return 0
    require_distinct && (link1 == link3 || link2 == link4) && return 0

    link3 > 0 && _shull_replace_neighbour!(triads[link3], neighbour_index, index)
    link2 > 0 && _shull_replace_neighbour!(triads[link2], index, neighbour_index)
    triads[neighbour_index] = second
    triads[index] = first
    return neighbour_index
end

function _shull_flip_one!(points::Vector{_SHullPoint}, triads::Vector{_SHullTriad},
                          index::Int; boundary_only::Bool=false,
                          require_distinct::Bool=true)
    triad = triads[index]
    for edge in 1:3
        neighbour = edge == 1 ? triad.bc : edge == 2 ? triad.ab : triad.ac
        neighbour == 0 && continue
        if boundary_only
            on_boundary = edge == 1 ? (triad.ac == 0 || triad.ab == 0) :
                          edge == 2 ? (triad.ac == 0 || triad.bc == 0) :
                                      (triad.bc == 0 || triad.ab == 0)
            on_boundary || continue
        end
        flipped = _shull_try_flip!(points, triads, index, edge; require_distinct)
        flipped > 0 && return flipped
    end
    return 0
end

function _shull_flip_pass!(points::Vector{_SHullPoint}, triads::Vector{_SHullTriad},
                           indices; boundary_only::Bool=false,
                           require_distinct::Bool=true)
    changed = Int[]
    for index in indices
        neighbour = _shull_flip_one!(points, triads, index; boundary_only, require_distinct)
        neighbour > 0 && append!(changed, (index, neighbour))
    end
    return changed
end

"""Complete SHull's three source-order Cline-Renka flip phases in place.

This directly follows `T_flip_pro`, `T_flip_pro_idx`, and `T_flip_edge` from
`interp`'s SHull source. It remains an internal prototype until its ordered
triads and downstream alpha-hull arcs agree with the R reference.
"""
function _shull_flip_triads!(state)
    points_by_id = Vector{_SHullPoint}(undef, length(state.points))
    for point in state.points
        points_by_id[point.id] = point
    end
    triads = state.triads
    pending = _shull_flip_pass!(points_by_id, triads, eachindex(triads))
    iteration = 1
    maximum_pending = 2 * (2 * length(state.points) - length(state.hull) - 2)
    while !isempty(pending) && iteration < 50
        pending = _shull_flip_pass!(points_by_id, triads, pending)
        length(pending) > maximum_pending && throw(ErrorException("SHull flip queue grew unexpectedly"))
        iteration += 1
    end

    pending = _shull_flip_pass!(
        points_by_id,
        triads,
        eachindex(triads);
        boundary_only=true,
        require_distinct=false,
    )
    iteration = 0
    while !isempty(pending) && iteration < 100
        pending = _shull_flip_pass!(points_by_id, triads, pending)
        iteration += 1
    end
    return state
end

"""Return final SHull triads after deterministic Float32 flip resolution.

This experimental function is intentionally separate from `delvor` while the
R-compatible triad and arc ordering is validated.
"""
function _shull_final_triads(xy::AbstractMatrix; rng=Random.MersenneTwister(0))
    try
        state = _shull_incremental_triads(xy)
        return _shull_flip_triads!(state)
    catch error
        _shull_requires_jitter(error) || rethrow()
        # R retries the topology only, then restores its original x/y vectors
        # before exposing tri.mesh. Its C++ trlist is nevertheless oriented
        # against these Float64 jitter coordinates, so retain them separately.
        jittered_xy = _shull_jittered_points(xy; rng=rng)
        state = _shull_incremental_triads(jittered_xy)
        return (; _shull_flip_triads!(state)..., orientation_xy=jittered_xy)
    end
end

"""Build R-compatible `interp::tri.mesh()\$trlist` rows from final SHull triads.

Julia's internal triads require the established orientation conversion to
match the C++ output. After a jitter retry, `tri.mesh()` restores public x/y
but C++ has already oriented `trlist` with its Float64 jitter input. The
working SHull points are Float32 and cannot reproduce this distinction, so
the retained `orientation_xy` is used only for this export step. Arc indices
are assigned in the same first-seen order as `shullDeltri.cpp`.
"""
function _shull_trlist(xy::AbstractMatrix, state=_shull_final_triads(xy))
    size(xy, 2) == 2 || throw(ArgumentError("xy must have two columns"))
    rows = Matrix{Int}(undef, length(state.triads), 9)
    orientation_xy = hasproperty(state, :orientation_xy) ? state.orientation_xy : xy
    for (row, triad) in enumerate(state.triads)
        orientation = (orientation_xy[triad.c, 1] - orientation_xy[triad.b, 1]) *
                      (orientation_xy[triad.b, 2] - orientation_xy[triad.a, 2]) +
                      (orientation_xy[triad.c, 2] - orientation_xy[triad.b, 2]) *
                      (orientation_xy[triad.a, 1] - orientation_xy[triad.b, 1])
        if orientation < 0
            rows[row, 1:6] .= (triad.a, triad.b, triad.c, triad.bc, triad.ac, triad.ab)
        else
            rows[row, 1:6] .= (triad.a, triad.c, triad.b, triad.bc, triad.ab, triad.ac)
        end
    end
    arc_indices = Dict{Tuple{Int, Int}, Int}()
    next_arc = 1
    for row in axes(rows, 1)
        for (column, left, right) in (
            (7, rows[row, 2], rows[row, 3]),
            (8, rows[row, 3], rows[row, 1]),
            (9, rows[row, 1], rows[row, 2]),
        )
            key = left < right ? (left, right) : (right, left)
            arc = get(arc_indices, key, 0)
            if arc == 0
                arc_indices[key] = next_arc
                arc = next_arc
                next_arc += 1
            end
            rows[row, column] = arc
        end
    end
    return rows
end

"""Reproduce the centre calculation used by `interp::circum()`.

The side lengths and barycentric reconstruction are Float64, while Heron's
semiperimeter, area, and radii are C++ `float`. Although only the centre is
needed by `alphahull::delvor()`, preserving those intermediate conversions
keeps its degenerate-case behaviour available to the SHull adapter.
"""
function _shull_interp_circumcentre(x1::Float64, y1::Float64,
                                    x2::Float64, y2::Float64,
                                    x3::Float64, y3::Float64)
    a = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    b = sqrt((x3 - x2)^2 + (y3 - y2)^2)
    c = sqrt((x1 - x3)^2 + (y1 - y3)^2)
    semi_perimeter = Float32((a + b + c) / 2.0)
    heron = semi_perimeter * (semi_perimeter - a) *
            (semi_perimeter - b) * (semi_perimeter - c)
    # C++ std::sqrt returns NaN for a slightly negative near-collinear Heron
    # term. Julia throws on a negative real input, so make that C++ behaviour
    # explicit. The centre below uses the Float64 side lengths directly.
    area = Float32(heron < 0 ? NaN : sqrt(heron))
    radius = Float32(a * b * c / (4.0 * area))
    inradius = area / semi_perimeter
    aspect_ratio = inradius / radius

    p = a^2 * (-a^2 + b^2 + c^2)
    q = b^2 * (a^2 - b^2 + c^2)
    r = c^2 * (a^2 + b^2 - c^2)
    total = p + q + r
    p /= total
    q /= total
    r /= total
    return q * x1 + r * x2 + p * x3, q * y1 + r * y2 + p * y3, aspect_ratio
end

"""Reconstruct the exposed interp SHull boundary-chain vector.

R dummycoor calls in.convex.hull on the public triSht object, not SHull's
internal hull state. For near-degenerate data those can differ because the R
boundary-chain walk retains repeated vertices. Preserve that source behavior
because it selects the side of each dummy circumcentre.
"""
function _shull_boundary_chain(trlist::AbstractMatrix{<:Integer})
    starts = Int[]
    ends = Int[]
    for row in axes(trlist, 1)
        trlist[row, 4] == 0 && (push!(starts, trlist[row, 2]); push!(ends, trlist[row, 3]))
        trlist[row, 5] == 0 && (push!(starts, trlist[row, 3]); push!(ends, trlist[row, 1]))
        trlist[row, 6] == 0 && (push!(starts, trlist[row, 1]); push!(ends, trlist[row, 2]))
    end
    count = length(starts)
    count == 0 && return Int[]

    chain = Vector{Int}(undef, count)
    chain[1], chain[2] = starts[1], ends[1]
    current = 2
    while current < count
        found = false
        for index in 2:count
            if chain[current] == starts[index]
                chain[current + 1] = ends[index]
                found = true
            end
        end
        if !found
            for index in 2:count
                if chain[current] == ends[index]
                    chain[current + 1] = starts[index]
                    found = true
                end
            end
        end
        current += 1
    end
    return chain
end

function _shull_in_convex_hull(points::AbstractMatrix, hull::AbstractVector{<:Integer},
                               x::Float64, y::Float64;
                               eps::Float64=1e-16)
    # Match R in.convex.hull against its reconstructed public boundary chain.
    count = length(hull)
    @inbounds for index in 1:count
        first = hull[index]
        second = hull[mod1(index + 1, count)]
        x1, y1 = points[first, 1], points[first, 2]
        x2, y2 = points[second, 1], points[second, 2]
        (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1) >= eps || return false
    end
    return true
end

function _shull_dummycoor(points::AbstractMatrix, first::Int, second::Int,
                          centre::NTuple{2,Float64}, away::Float64, hull)
    # This reproduces alphahull::dummycoor(). In particular, infer the
    # exterior side from R's global convex-hull test, not from the third
    # vertex of one (possibly near-degenerate) boundary triangle.
    vx = points[second, 2] - points[first, 2]
    vy = -(points[second, 1] - points[first, 1])
    norm2 = vx^2 + vy^2
    norm2 > 0 && ((vx, vy) = (vx / norm2, vy / norm2))
    midpoint_x = (points[first, 1] + points[second, 1]) / 2
    midpoint_y = (points[first, 2] + points[second, 2]) / 2
    in_hull = _shull_in_convex_hull(
        points,
        hull,
        midpoint_x + 1e-5 * vx,
        midpoint_y + 1e-5 * vy,
    )
    return in_hull ? (centre[1] - away * vx, centre[2] - away * vy) :
                     (centre[1] + away * vx, centre[2] + away * vy)
end

"""Build the R-compatible twelve-column `alphahull::delvor()` mesh.

The result uses the existing Julia `DelVor` mesh schema. It is intentionally
kept separate from `delvor` until differential tests cover degenerate inputs
and downstream alpha-hull construction.
"""
function _shull_mesh(xy::AbstractMatrix, state=_shull_final_triads(xy))
    points = Matrix{Float64}(xy)
    trlist = _shull_trlist(points, state)
    boundary_chain = _shull_boundary_chain(trlist)
    triangle_count = size(trlist, 1)
    centres = Vector{NTuple{2, Float64}}(undef, triangle_count)
    for row in axes(trlist, 1)
        a, b, c = trlist[row, 1], trlist[row, 2], trlist[row, 3]
        centre_x, centre_y, _ = _shull_interp_circumcentre(
            points[a, 1], points[a, 2], points[b, 1], points[b, 2], points[c, 1], points[c, 2],
        )
        centres[row] = (centre_x, centre_y)
    end

    # This is the exact aux1/aux2/aux3 rbind order in alphahull::delvor().
    auxiliary = NTuple{5, Int}[]
    for row in axes(trlist, 1)
        push!(auxiliary, (trlist[row, 7], trlist[row, 2], trlist[row, 3], row, trlist[row, 4]))
    end
    for row in axes(trlist, 1)
        push!(auxiliary, (trlist[row, 8], trlist[row, 1], trlist[row, 3], row, trlist[row, 5]))
    end
    for row in axes(trlist, 1)
        push!(auxiliary, (trlist[row, 9], trlist[row, 1], trlist[row, 2], row, trlist[row, 6]))
    end
    seen_arcs = Set{Int}()
    filter!(row -> row[1] ∉ seen_arcs && (push!(seen_arcs, row[1]); true), auxiliary)

    away = max(maximum(points[:, 1]) - minimum(points[:, 1]),
               maximum(points[:, 2]) - minimum(points[:, 2]))
    mesh = Matrix{Float64}(undef, length(auxiliary), 12)
    for (row, (arc, first, second, first_triangle, second_triangle)) in enumerate(auxiliary)
        first_centre = centres[first_triangle]
        second_centre = if second_triangle == 0
            dummy = _shull_dummycoor(
                points,
                first,
                second,
                first_centre,
                away,
                boundary_chain,
            )
            push!(centres, dummy)
            dummy
        else
            centres[second_triangle]
        end
        mesh[row, :] .= (
            first,
            second,
            points[first, 1],
            points[first, 2],
            points[second, 1],
            points[second, 2],
            first_centre[1],
            first_centre[2],
            second_centre[1],
            second_centre[2],
            first_triangle == 0,
            second_triangle == 0,
        )
    end
    return mesh
end
