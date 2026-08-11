"""Construct the beta-skeleton style disk-window arcs used by `alphahull::dw`.

The returned matrix has the same nine columns as the R implementation:
`c1, c2, r, v.x, v.y, theta, point, orig.point, inte`.
Indices are one-based, as in the rest of this Julia package.
"""
function dw(x, y=nothing; eps)
    eps > 0 || throw(ArgumentError("eps must be positive"))
    sample = _points(x, y)
    extremes = ashape(sample; alpha=eps).alpha_extremes
    n = size(sample, 1)
    m = length(extremes)
    d = zeros(Float64, m, n)
    vx = zeros(Float64, m, n)
    vy = zeros(Float64, m, n)
    theta = zeros(Float64, m, n)
    for (row, idx) in enumerate(extremes)
        vx[row, :] .= sample[:, 1] .- sample[idx, 1]
        vy[row, :] .= sample[:, 2] .- sample[idx, 2]
        d[row, :] .= hypot.(vx[row, :], vy[row, :])
        nz = d[row, :] .> 0
        vx[row, nz] ./= d[row, nz]
        vy[row, nz] ./= d[row, nz]
        hit = nz .& (d[row, :] .<= 2eps)
        theta[row, hit] .= acos.(clamp.(d[row, hit] .* (0.5 / eps), -1.0, 1.0))
    end

    rows = Vector{Vector{Float64}}()
    for (row, idx) in enumerate(extremes)
        push!(rows, [sample[idx, 1], sample[idx, 2], eps, 1.0, 0.0, pi/2,
                     Float64(row), Float64(idx), 0.0])
    end
    for (row, idx) in enumerate(extremes)
        push!(rows, [sample[idx, 1], sample[idx, 2], eps, -1.0, 0.0, pi/2,
                     Float64(row), Float64(idx), 0.0])
    end
    nowatch = Int[]
    watch = 1
    while watch <= length(rows)
        row = rows[watch]
        prow = Int(row[7])
        candidates = findall(j -> d[prow, j] > 0 && d[prow, j] <= 2eps, 1:n)
        for jidx in candidates
            avx, avy = row[4], row[5]
            angox = avy >= 0 ? acos(clamp(avx, -1.0, 1.0)) : 2pi - acos(clamp(avx, -1.0, 1.0))
            vintx, vinty = vx[prow, jidx], vy[prow, jidx]
            rvx, rvy = rotation([vintx, vinty], angox)
            intertheta = theta[prow, jidx]
            if rvy >= 0
                ang = acos(clamp(rvx, -1.0, 1.0))
                vals = [-row[6], row[6], ang-intertheta, ang+intertheta]
                labels = (:theta1, :theta2, :beta1, :beta2)
                order = Tuple(labels[sortperm(vals)])
                if order == (:theta1, :beta1, :theta2, :beta2)
                    middle = (vals[3] - vals[1]) / 2
                    nvx, nvy = rotation([1.0, 0.0], row[6] - middle - angox)
                    row[4:6] = [nvx, nvy, middle]; row[9] = jidx
                elseif order == (:beta1, :theta1, :theta2, :beta2)
                    push!(nowatch, watch)
                elseif order == (:theta1, :beta1, :beta2, :theta2)
                    middle = (vals[3] - vals[1]) / 2
                    nvx, nvy = rotation([1.0, 0.0], row[6] - middle - angox)
                    middle2 = (vals[2] - vals[4]) / 2
                    nvx2, nvy2 = rotation([1.0, 0.0], -row[6] + middle2 - angox)
                    push!(rows, [row[1], row[2], row[3], nvx2, nvy2, middle2,
                                 row[7], row[8], Float64(jidx)])
                    row[4:6] = [nvx, nvy, middle]; row[9] = jidx
                elseif order == (:beta1, :theta1, :beta2, :theta2)
                    middle2 = (vals[2] - vals[4]) / 2
                    nvx2, nvy2 = rotation([1.0, 0.0], -row[6] + middle2 - angox)
                    row[4:6] = [nvx2, nvy2, middle2]; row[9] = jidx
                end
            else
                ang = acos(clamp(rvx, -1.0, 1.0))
                vals = [-row[6], row[6], -ang-intertheta, -ang+intertheta]
                labels = (:theta1, :theta2, :beta1, :beta2)
                order = Tuple(labels[sortperm(vals)])
                if order == (:theta1, :beta1, :theta2, :beta2)
                    middle = (vals[3] - vals[1]) / 2
                    nvx, nvy = rotation([1.0, 0.0], row[6] - middle - angox)
                    row[4:6] = [nvx, nvy, middle]; row[9] = jidx
                elseif order == (:beta1, :theta1, :theta2, :beta2)
                    push!(nowatch, watch)
                elseif order == (:theta1, :beta1, :beta2, :theta2)
                    middle = (vals[3] - vals[1]) / 2
                    nvx, nvy = rotation([1.0, 0.0], row[6] - middle - angox)
                    middle2 = (vals[2] - vals[4]) / 2
                    nvx2, nvy2 = rotation([1.0, 0.0], -row[6] + middle2 - angox)
                    push!(rows, [row[1], row[2], row[3], nvx2, nvy2, middle2,
                                 row[7], row[8], Float64(jidx)])
                    row[4:6] = [nvx, nvy, middle]; row[9] = jidx
                elseif order == (:beta1, :theta1, :beta2, :theta2)
                    middle2 = (vals[2] - vals[4]) / 2
                    nvx2, nvy2 = rotation([1.0, 0.0], -row[6] + middle2 - angox)
                    row[4:6] = [nvx2, nvy2, middle2]; row[9] = jidx
                end
            end
        end
        watch += 1
    end
    for idx in sort(unique(nowatch); rev=true)
        1 <= idx <= length(rows) && deleteat!(rows, idx)
    end
    return isempty(rows) ? zeros(Float64, 0, 9) : reduce(vcat, permutedims.(rows))
end

"""Sample an arc as an `n`-by-2 matrix of coordinates."""
function arc(c, r::Real, v, theta::Real; n::Integer=100)
    n >= 2 || throw(ArgumentError("n must be at least two"))
    angles = anglesArc(v, theta)
    t = range(angles[1], angles[2]; length=n)
    return hcat(c[1] .+ r .* cos.(t), c[2] .+ r .* sin.(t))
end

function _track_points(x, y, nps, sc, rng)
    sample = _points(x, y)
    n = size(sample, 1)
    nps >= 0 || throw(ArgumentError("nps must be nonnegative"))
    lengths = [hypot(sample[i+1,1]-sample[i,1], sample[i+1,2]-sample[i,2]) for i in 1:n-1]
    total = sum(lengths)
    total > 0 || throw(ArgumentError("track must contain at least one nonzero segment"))
    extra = Matrix{Float64}(undef, nps, 2)
    cumulative_lengths = cumsum(lengths)
    if nps > 0
        for k in 1:nps
            z = rand(rng) * total
            seg = findfirst(cumulative_lengths .>= z)
            seg = isnothing(seg) ? n-1 : seg
            u = rand(rng)
            extra[k, :] .= (sample[seg, :] .* (1-u) .+ sample[seg+1, :] .* u)
        end
    end
    allpts = vcat(sample, extra)
    jitter = max(eps(Float64), 1e-10 * max(1.0, maximum(abs, allpts)))
    return (allpts .+ jitter .* randn(rng, size(allpts))) .* sc
end

"""Generate sampled disk-window arc paths for a track; returns a vector of matrices.

Unlike R's plotting helper, this returns the sampled coordinates directly.
"""
function dw_track(x, y=nothing; eps, nps::Integer=20000, sc::Real=100, rng=Random.default_rng())
    p = _track_points(x, y, nps, sc, rng)
    w = dw(p; eps=eps*sc)
    return [arc(w[i,1:2], w[i,3], w[i,4:5], w[i,6]; n=100) ./ sc for i in axes(w,1)]
end

"""Generate sampled alpha-hull arc paths for a track; returns a vector of matrices.

Unlike R's plotting helper, this returns the sampled coordinates directly.
"""
function ahull_track(x, y=nothing; alpha, nps::Integer=10000, sc::Real=100, rng=Random.default_rng())
    p = _track_points(x, y, nps, sc, rng)
    h = ahull(p; alpha=alpha*sc)
    rows = findall(h.arcs[:,3] .> 0)
    return [arc(h.arcs[i,1:2], h.arcs[i,3], h.arcs[i,4:5], h.arcs[i,6]; n=100) ./ sc for i in rows]
end
