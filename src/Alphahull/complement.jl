"""Return the 19-column complement table used by R `alphahull::complement`."""
function complement(x, y=nothing; alpha)
    alpha >= 0 || throw(ArgumentError("alpha must be nonnegative"))
    dv = x isa DelVor && y === nothing ? x : delvor(x, y)
    m = dv.mesh
    # `alphahull::complement()` derives `vert <- mesh[, "mx1"] ==
    # mesh[, "mx2"]` before constructing any balls.  An NA circumcentre
    # makes `sum(vert)` NA and R aborts the whole candidate.  SHull exposes
    # that same C++ state as `NaN`; do not let Julia's boolean comparisons
    # silently turn it into a usable alpha hull.
    any(isnan, @view m[:, 7:10]) && throw(ArgumentError(
        "alpha-hull complement is undefined for a NaN circumcentre",
    ))
    ne = size(m, 1)
    dm1 = hypot.(m[:,3] .- m[:,7], m[:,4] .- m[:,8])
    dm2 = hypot.(m[:,3] .- m[:,9], m[:,4] .- m[:,10])
    pmx, pmy = (m[:,3] .+ m[:,5]) ./ 2, (m[:,4] .+ m[:,6]) ./ 2
    dm = hypot.(m[:,3] .- m[:,5], m[:,4] .- m[:,6]) ./ 2
    d1 = hypot.(m[:,7] .- pmx, m[:,8] .- pmy)
    d2 = hypot.(m[:,9] .- pmx, m[:,10] .- pmy)
    theta1 = atan.(dm ./ d1)
    theta2 = atan.(dm ./ d2)
    vx1 = ifelse.(d1 .!= 0, (pmx .- m[:,7]) ./ d1, 0.0)
    vy1 = ifelse.(d1 .!= 0, (pmy .- m[:,8]) ./ d1, 0.0)
    vx2 = ifelse.(d2 .!= 0, (pmx .- m[:,9]) ./ d2, 0.0)
    vy2 = ifelse.(d2 .!= 0, (pmy .- m[:,10]) ./ d2, 0.0)
    betw = falses(ne)
    for i in 1:ne
        betw[i] = m[i,7] == m[i,9] ? min(m[i,8],m[i,10]) <= pmy[i] <= max(m[i,8],m[i,10]) :
                                      min(m[i,7],m[i,9]) <= pmx[i] <= max(m[i,7],m[i,9])
    end
    aux = alpha^2 .- dm.^2
    rows = NTuple{19,Float64}[]
    function pushball!(i, c1, c2, r, ind, vx, vy, theta)
        push!(rows, (c1, c2, r,
                     m[i, 1], m[i, 2], m[i, 3], m[i, 4], m[i, 5], m[i, 6],
                     m[i, 7], m[i, 8], m[i, 9], m[i, 10], m[i, 11], m[i, 12],
                     Float64(ind), vx, vy, theta))
    end
    for i in 1:ne
        if m[i,11] == 0 && dm1[i] > alpha
            pushball!(i, m[i,7], m[i,8], dm1[i], 1, vx1[i], vy1[i], theta1[i])
        end
        a = aux[i] > 0 ? sqrt(aux[i]) : 0.0
        if aux[i] > 0 && betw[i] && (m[i,11] == 0 ? dm1[i] > alpha : true)
            pushball!(i, pmx[i]-a*vx1[i], pmy[i]-a*vy1[i], alpha, 1, vx1[i], vy1[i], atan(dm[i]/a))
        elseif aux[i] > 0 && !betw[i] &&
               (m[i,11] == 0 ? dm1[i] > alpha : true) && alpha-a < dm2[i]-d2[i]
            pushball!(i, pmx[i]-a*vx1[i], pmy[i]-a*vy1[i], alpha, 1, vx1[i], vy1[i], atan(dm[i]/a))
        end
        if m[i,12] == 0 && dm2[i] > alpha
            pushball!(i, m[i,9], m[i,10], dm2[i], 2, vx2[i], vy2[i], theta2[i])
        end
        a = aux[i] > 0 ? sqrt(aux[i]) : 0.0
        if aux[i] > 0 && betw[i] && (m[i,12] == 0 ? dm2[i] > alpha : true)
            pushball!(i, pmx[i]-a*vx2[i], pmy[i]-a*vy2[i], alpha, 2, vx2[i], vy2[i], atan(dm[i]/a))
        elseif aux[i] > 0 && !betw[i] &&
               (m[i,12] == 0 ? dm2[i] > alpha : true) && alpha-a < dm1[i]-d1[i]
            pushball!(i, pmx[i]-a*vx2[i], pmy[i]-a*vy2[i], alpha, 2, vx2[i], vy2[i], atan(dm[i]/a))
        end
        if m[i,12] == 1
            if m[i,3] == m[i,5]
                a = m[i,3]
                sig = m[i,9] <= a ? -4.0 : -3.0
                pushball!(i, a, 0.0, sig, 2, 0.0, 0.0, 0.0)
            else
                b = (m[i,6]-m[i,4])/(m[i,5]-m[i,3])
                a = m[i,4]-m[i,3]*b
                sig = m[i,10] <= a+b*m[i,9] ? -2.0 : -1.0
                pushball!(i, a, b, sig, 2, 0.0, 0.0, 0.0)
            end
        end
    end
    isempty(rows) && return zeros(Float64, 0, 19)
    out = Matrix{Float64}(undef, length(rows), 19)
    for i in eachindex(rows)
        values = rows[i]
        @inbounds for j in 1:19
            out[i, j] = values[j]
        end
    end
    return out
end
