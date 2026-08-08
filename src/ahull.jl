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

function _angle_data(vx, vy, theta, ivx, ivy, itheta)
    angox = vy >= 0 ? acos(clamp(vx, -1.0, 1.0)) : 2pi - acos(clamp(vx, -1.0, 1.0))
    _, iy = _rotate(ivx, ivy, angox)
    ix, _ = _rotate(ivx, ivy, angox)
    phi = acos(clamp(ix, -1.0, 1.0))
    signedphi = iy >= 0 ? phi : -phi
    angles = [-theta, theta, signedphi - itheta, signedphi + itheta]
    labels = (:theta1, :theta2, :beta1, :beta2)
    order = Tuple(labels[sortperm(angles)])
    return (angles=angles, order=order, angox=angox,
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

function _initial_arcs(sh, cp)
    arcs = Vector{Vector{Float64}}()
    ends = NTuple{2,Int}[]
    for e in axes(sh.edges, 1)
        i, j = Int(sh.edges[e, 1]), Int(sh.edges[e, 2])
        matches = findall(k -> Int(cp[k, 4]) == i && Int(cp[k, 5]) == j, axes(cp, 1))
        isempty(matches) && continue
        first_side = count(k -> cp[k, 16] == 1, matches)
        0 < first_side < length(matches) && continue
        balls = filter(k -> cp[k, 3] > 0, matches)
        isempty(balls) && continue
        k = balls[argmin(cp[balls, 3])]
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
    watch = 1
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
                        case = ad.order == (:theta1, :beta1, :beta2, :theta2) ? 2 :
                               ad.order == (:beta1, :beta2, :theta1, :theta2) ? 1 : 0
                        if ad.order == (:beta1, :theta1, :beta2, :theta2)
                            case = abs((ad.angles[1] - ad.angles[3]) / 2) < 1e-5 ? 2 : 1
                        end
                        if case == 2
                            middle = (ad.angles[2] - ad.angles[4]) / 2
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
                        case = ad.order == (:theta1, :beta1, :beta2, :theta2) ? 2 :
                               ad.order == (:theta1, :theta2, :beta1, :beta2) ? 1 : 0
                        if ad.order == (:theta1, :beta1, :theta2, :beta2)
                            case = abs((ad.angles[2] - ad.angles[4]) / 2) < 1e-5 ? 2 : 1
                        end
                        if case == 2
                            middle = (ad.angles[3] - ad.angles[1]) / 2
                            nvx, nvy = _rotate(1.0, 0.0, aw[6] - middle - ad.angox)
                            qx, qy = _rotate(nvx, nvy, -middle)
                            push!(points, (aw[1] + aw[3] * qx, aw[2] + aw[3] * qy))
                            inn = length(points)
                            aw[4], aw[5], aw[6] = nvx, nvy, middle
                            ends[watch] = (ends[watch][1], inn)
                            _arc_from_endpoints!(aj, points, inn, ends[j][2])
                            ends[j] = (inn, ends[j][2])
                        end
                    elseif !shared && ad.order == (:theta1, :beta1, :beta2, :theta2)
                        bd = _angle_data(aj[4], aj[5], aj[6], ci.v2[1], ci.v2[2], ci.theta2)
                        if bd.order == (:theta1, :beta1, :beta2, :theta2)
                            middle1 = (ad.angles[3] - ad.angles[1]) / 2
                            nvx1, nvy1 = _rotate(1.0, 0.0, aw[6] - middle1 - ad.angox)
                            middle2 = (ad.angles[2] - ad.angles[4]) / 2
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
            end
            j += 1
        end
        watch += 1
    end
    order = Int[]
    remaining = collect(1:length(ends))
    while !isempty(remaining)
        current = popfirst!(remaining)
        push!(order, current)
        while true
            nextpos = findfirst(k -> ends[k][1] == ends[current][2], remaining)
            isnothing(nextpos) && break
            current = remaining[nextpos]
            deleteat!(remaining, nextpos)
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
