struct AHull
    arcs::Matrix{Float64}
    xahull::Matrix{Float64}
    length::Float64
    complement::Matrix{Float64}
    alpha::Float64
    ashape::AShape
end

_rotate(vx, vy, theta) = (cos(theta) * vx + sin(theta) * vy, -sin(theta) * vx + cos(theta) * vy)
lengthahull(arcs::AbstractMatrix) = sum(2 .* arcs[:, 6] .* arcs[:, 3])

const _ORDER_1324 = UInt16(1 | (3 << 3) | (4 << 6) | (2 << 9))
const _ORDER_3412 = UInt16(3 | (4 << 3) | (1 << 6) | (2 << 9))
const _ORDER_3142 = UInt16(3 | (1 << 3) | (4 << 6) | (2 << 9))
const _ORDER_1234 = UInt16(1 | (2 << 3) | (3 << 6) | (4 << 9))
const _ORDER_1324_ALT = UInt16(1 | (3 << 3) | (2 << 6) | (4 << 9))

@inline function _angle_order(a1, a2, a3, a4)
    v1, v2, v3, v4 = a1, a2, a3, a4
    i1, i2, i3, i4 = 1, 2, 3, 4
    if v2 < v1
        v1, v2 = v2, v1
        i1, i2 = i2, i1
    end
    if v4 < v3
        v3, v4 = v4, v3
        i3, i4 = i4, i3
    end
    if v3 < v1
        v1, v3 = v3, v1
        i1, i3 = i3, i1
    end
    if v4 < v2
        v2, v4 = v4, v2
        i2, i4 = i4, i2
    end
    if v3 < v2
        i2, i3 = i3, i2
    end
    return UInt16(i1 | (i2 << 3) | (i3 << 6) | (i4 << 9))
end

@inline function _angle_data(vx, vy, theta, ivx, ivy, itheta)
    angox = vy >= 0 ? acos(clamp(vx, -1.0, 1.0)) : 2pi - acos(clamp(vx, -1.0, 1.0))
    _, iy = _rotate(ivx, ivy, angox)
    ix, _ = _rotate(ivx, ivy, angox)
    phi = acos(clamp(ix, -1.0, 1.0))
    signedphi = iy >= 0 ? phi : -phi
    angle1, angle2 = -theta, theta
    angle3, angle4 = signedphi - itheta, signedphi + itheta
    return (angle1=angle1, angle2=angle2, angle3=angle3, angle4=angle4,
            order=_angle_order(angle1, angle2, angle3, angle4), angox=angox,
            theta1=angox-theta, theta2=angox+theta,
            beta1=signedphi+angox-itheta, beta2=signedphi+angox+itheta)
end

function _arc_from_endpoints!(arc, points, a, b)
    px, py = points[a]
    qx, qy = points[b]
    mx, my = (px + qx) / 2, (py + qy) / 2
    dx, dy = mx - arc[1], my - arc[2]
    d = hypot(dx, dy)
    arc[4], arc[5] = dx / d, dy / d
    arc[6] = atan(hypot(mx - px, my - py) / d)
end

mutable struct _ComplementArcSummary
    count::Int
    first_side_count::Int
    smallest_ball::Int
end

function _initial_arcs(sh, cp)
    summaries = Dict{Tuple{Int,Int},_ComplementArcSummary}()
    for k in axes(cp, 1)
        key = (Int(cp[k, 4]), Int(cp[k, 5]))
        summary = get!(summaries, key) do
            _ComplementArcSummary(0, 0, 0)
        end
        summary.count += 1
        summary.first_side_count += cp[k, 16] == 1
        if cp[k, 3] > 0 &&
           (summary.smallest_ball == 0 || cp[k, 3] < cp[summary.smallest_ball, 3])
            summary.smallest_ball = k
        end
    end

    arcs = Vector{Vector{Float64}}()
    ends = NTuple{2,Int}[]
    for e in axes(sh.edges, 1)
        i, j = Int(sh.edges[e, 1]), Int(sh.edges[e, 2])
        summary = get(summaries, (i, j), nothing)
        isnothing(summary) && continue
        0 < summary.first_side_count < summary.count && continue
        k = summary.smallest_ball
        # R `ahull()` calls `min()` on the positive-radius complement rows
        # for every selected alpha-shape edge. With no such row it errors
        # (`replacement has length zero`) rather than silently omitting the
        # edge. That failure is part of rangeBuilder's alpha retry path.
        k == 0 && throw(ArgumentError(
            "alpha-hull complement has no positive radius for an R-selected edge",
        ))
        arc = [cp[k, 1], cp[k, 2], cp[k, 3], cp[k, 17], cp[k, 18], cp[k, 19]]
        pmx, pmy = (cp[k, 6] + cp[k, 8]) / 2, (cp[k, 7] + cp[k, 9]) / 2
        tx, ty = _rotate(arc[4], arc[5], arc[6])
        a2 = (cp[k, 6] - pmx) * tx + (cp[k, 7] - pmy) * ty
        push!(arcs, arc)
        push!(ends, a2 > 0 ? (i, j) : (j, i))
    end
    return arcs, ends
end

"""Port of the arc cutting and ordering state machine in R `ahull`.

The mutable `arcs`, `ends`, and `points` collections correspond respectively
to the R variables `arcs`, `indp`, and `cutp`.
"""
function _cut_and_order!(arcs, ends, points)
    # R walks a dynamically growing watch-by-j matrix. Preserving its source
    # order is necessary because each cut can change later endpoint states.
    watch = 1
    case = 0
    case_defined = false
    while watch <= length(arcs)
        j = 1
        while j <= length(arcs)
            if j != watch
                aw, aj = arcs[watch], arcs[j]
                ci = inter(aw[1], aw[2], aw[3], aj[1], aj[2], aj[3])
                if ci.n_cut == 2
                    ad = _angle_data(aw[4], aw[5], aw[6], ci.v1[1], ci.v1[2], ci.theta1)
                    shared = ends[watch][1] == ends[j][1] || ends[watch][1] == ends[j][2] ||
                             ends[watch][2] == ends[j][1] || ends[watch][2] == ends[j][2]
                    if shared && ends[watch][1] == ends[j][2]
                        next_case = ad.order == _ORDER_1324 ? 2 :
                                    ad.order == _ORDER_3412 ? 1 : 0
                        if ad.order == _ORDER_3142
                            next_case = abs((ad.angle1 - ad.angle3) / 2) < 1e-5 ? 2 : 1
                        end
                        # R assigns `case <- 0` only after an earlier
                        # `j != watch` iteration. An unsupported first pair
                        # errors; later unsupported pairs retain zero.
                        next_case != 0 && (case = next_case)
                        case == 0 && !case_defined && throw(ArgumentError(
                            "alpha-hull arc cut has an undefined R case",
                        ))
                        if case == 2
                            middle = (ad.angle2 - ad.angle4) / 2
                            nvx, nvy = _rotate(1.0, 0.0, -aw[6] + middle - ad.angox)
                            qx, qy = _rotate(nvx, nvy, middle)
                            push!(points, (aw[1] + aw[3] * qx, aw[2] + aw[3] * qy))
                            inn = length(points)
                            aw[4], aw[5], aw[6] = nvx, nvy, middle
                            ends[watch] = (inn, ends[watch][2])
                            _arc_from_endpoints!(aj, points, inn, ends[j][1])
                            ends[j] = (ends[j][1], inn)
                        end
                    elseif shared && ends[watch][2] == ends[j][1]
                        next_case = ad.order == _ORDER_1324 ? 2 :
                                    ad.order == _ORDER_1234 ? 1 : 0
                        if ad.order == _ORDER_1324_ALT
                            next_case = abs((ad.angle2 - ad.angle4) / 2) < 1e-5 ? 2 : 1
                        end
                        next_case != 0 && (case = next_case)
                        case == 0 && !case_defined && throw(ArgumentError(
                            "alpha-hull arc cut has an undefined R case",
                        ))
                        if case == 2
                            middle = (ad.angle3 - ad.angle1) / 2
                            nvx, nvy = _rotate(1.0, 0.0, aw[6] - middle - ad.angox)
                            qx, qy = _rotate(nvx, nvy, -middle)
                            push!(points, (aw[1] + aw[3] * qx, aw[2] + aw[3] * qy))
                            inn = length(points)
                            aw[4], aw[5], aw[6] = nvx, nvy, middle
                            ends[watch] = (ends[watch][1], inn)
                            _arc_from_endpoints!(aj, points, inn, ends[j][2])
                            ends[j] = (inn, ends[j][2])
                        end
                    elseif !shared && ad.order == _ORDER_1324
                        bd = _angle_data(aj[4], aj[5], aj[6], ci.v2[1], ci.v2[2], ci.theta2)
                        if bd.order == _ORDER_1324
                            middle1 = (ad.angle3 - ad.angle1) / 2
                            nvx1, nvy1 = _rotate(1.0, 0.0, aw[6] - middle1 - ad.angox)
                            middle2 = (ad.angle2 - ad.angle4) / 2
                            nvx2, nvy2 = _rotate(1.0, 0.0, -aw[6] + middle2 - ad.angox)
                            push!(arcs, [aw[1], aw[2], aw[3], nvx2, nvy2, middle2])
                            aw[4], aw[5], aw[6] = nvx1, nvy1, middle1
                            q1x, q1y = _rotate(nvx1, nvy1, -middle1)
                            q2x, q2y = _rotate(nvx2, nvy2, middle2)
                            p1 = (aw[1] + aw[3] * q1x, aw[2] + aw[3] * q1y)
                            p2 = (aw[1] + aw[3] * q2x, aw[2] + aw[3] * q2y)
                            old_watch_end = ends[watch][2]
                            push!(points, p1); inn1 = length(points)
                            push!(points, p2); inn2 = length(points)
                            ends[watch] = (ends[watch][1], inn1)
                            push!(ends, (inn2, old_watch_end))
                            old_j_start = ends[j][1]
                            ends[j] = (inn1, ends[j][2])
                            push!(ends, (old_j_start, inn2))
                            _arc_from_endpoints!(aj, points, inn1, ends[j][2])
                            newarc = [aj[1], aj[2], aj[3], aj[4], aj[5], aj[6]]
                            _arc_from_endpoints!(newarc, points, inn2, old_j_start)
                            push!(arcs, newarc)
                        end
                    end
                    end
                case = 0
                case_defined = true
            end
            j += 1
        end
        watch += 1
    end
    order = Int[]
    sizehint!(order, length(ends))
    by_start = [Int[] for _ in 1:length(points)]
    for i in eachindex(ends)
        push!(by_start[ends[i][1]], i)
    end
    next_position = ones(Int, length(points))
    used = falses(length(ends))
    first_remaining = 1
    remaining_count = length(ends)
    while remaining_count > 0
        while used[first_remaining]
            first_remaining += 1
        end
        current = first_remaining
        used[current] = true
        remaining_count -= 1
        push!(order, current)
        while true
            endpoint = ends[current][2]
            candidates = by_start[endpoint]
            position = next_position[endpoint]
            while position <= length(candidates) && used[candidates[position]]
                position += 1
            end
            next_position[endpoint] = position
            position > length(candidates) && break
            current = candidates[position]
            next_position[endpoint] = position + 1
            used[current] = true
            remaining_count -= 1
            push!(order, current)
        end
    end
    return arcs[order], ends[order], points
end

"""Compute the alpha-convex hull with R-compatible arc cutting and ordering."""
function ahull(x, y=nothing; alpha)
    alpha >= 0 || throw(ArgumentError("alpha must be nonnegative"))
    dv = x isa DelVor && y === nothing ? x : delvor(x, y)
    sh = ashape(dv; alpha)
    cp = complement(dv; alpha)
    arcs, ends = _initial_arcs(sh, cp)
    points = [(sh.x[i, 1], sh.x[i, 2]) for i in axes(sh.x, 1)]
    if !isempty(arcs)
        arcs, ends, points = _cut_and_order!(arcs, ends, points)
    end
    used = Set{Int}()
    for endpoint in ends
        push!(used, endpoint[1]); push!(used, endpoint[2])
    end
    nextra = count(i -> !(i in used), sh.alpha_extremes)
    out = Matrix{Float64}(undef, length(arcs) + nextra, 8)
    for (row, (arc, (a, b))) in enumerate(zip(arcs, ends))
        @inbounds begin
            out[row, 1] = arc[1]
            out[row, 2] = arc[2]
            out[row, 3] = arc[3]
            out[row, 4] = arc[4]
            out[row, 5] = arc[5]
            out[row, 6] = arc[6]
            out[row, 7] = a
            out[row, 8] = b
        end
    end
    row = length(arcs)
    for i in sh.alpha_extremes
        if !(i in used)
            row += 1
            @inbounds begin
                out[row, 1] = sh.x[i, 1]
                out[row, 2] = sh.x[i, 2]
                out[row, 3] = 0.0
                out[row, 4] = 0.0
                out[row, 5] = 0.0
                out[row, 6] = 0.0
                out[row, 7] = i
                out[row, 8] = i
            end
        end
    end
    xahull = Matrix{Float64}(undef, length(points), 2)
    for (row, point) in enumerate(points)
        @inbounds begin
            xahull[row, 1] = point[1]
            xahull[row, 2] = point[2]
        end
    end
    return AHull(out, xahull, lengthahull(out), cp, Float64(alpha), sh)
end
