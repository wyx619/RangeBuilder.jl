function _point_in_polygon(p, poly)
    isempty(poly) && return false
    ring = Meshes.Ring([(poly[i, 1], poly[i, 2]) for i in axes(poly, 1)])
    return Meshes.sideof(Meshes.Point(p[1], p[2]), ring) != Meshes.OUT
end

_polygon_area(poly) = abs(sum(poly[i,1] * poly[i == size(poly,1) ? 1 : i+1,2] -
                              poly[i == size(poly,1) ? 1 : i+1,1] * poly[i,2]
                              for i in axes(poly, 1))) / 2

"""Evaluate the boundary-area construction used by R `areaahulleval`."""
function areaahulleval(h::AHull)
    arcs = copy(h.arcs)
    size(arcs, 1) == 0 && return 0.0
    # This is the literal pair-swapping pass used by areaahulleval.R.
    endpoint = Int.(arcs[:, 7:8])
    aa = vec(permutedims(endpoint))
    nuevoa = copy(aa)
    miropos, newpos = 2, 3
    while miropos < length(aa)
        positions = findall(==(aa[miropos]), aa)
        positions = filter(>(miropos), positions)
        isempty(positions) && break
        pos = first(positions)
        if iseven(pos)
            oldpair = (pos, pos - 1)
            filan = pos ÷ 2
        else
            oldpair = (pos, pos + 1)
            filan = pos ÷ 2 + 1
        end
        newpair = (newpos, newpos + 1)
        tmp = (nuevoa[newpair[1]], nuevoa[newpair[2]])
        nuevoa[newpair[1]], nuevoa[newpair[2]] = aa[oldpair[1]], aa[oldpair[2]]
        nuevoa[oldpair[1]], nuevoa[oldpair[2]] = tmp
        filold = newpos ÷ 2 + 1
        if filan != filold && filan >= 1 && filold <= size(arcs, 1)
            tmp_row = copy(arcs[filan, :])
            arcs[filan, :] = arcs[filold, :]
            arcs[filold, :] = tmp_row
        end
        aa = copy(nuevoa)
        miropos += 2
        newpos += 2
    end
    arcs[:, 7:8] = reshape(nuevoa, 2, :)'
    positive = findall(arcs[:, 3] .> 0)
    isempty(positive) && return 0.0
    arcs = arcs[positive, :]
    comps = zeros(Int, size(arcs, 1))
    row, component = 1, 0
    while row <= size(arcs, 1)
        check1 = arcs[row, 7]
        rownew = findfirst(==(check1), arcs[:, 8])
        isnothing(rownew) && break
        component += 1
        comps[row:rownew] .= component
        row = rownew + 1
    end
    # The ahull arc ordering is a collection of closed chains. Compute each
    # chain's polygonal shoelace term plus its circular-segment corrections.
    total = 0.0
    for comp in 1:maximum(comps)
        inds = findall(==(comp), comps)
        endpoint_ids = Int.(unique(vec(arcs[inds, 7:8])))
        poly = h.xahull[endpoint_ids, :]
        hole = -1.0
        small = findfirst(i -> arcs[inds[i], 6] < pi / 2, eachindex(inds))
        if !isnothing(small)
            c = arcs[inds[small], 1:2]
            hole = _point_in_polygon(c, poly) ? 1.0 : -1.0
        end
        circular = sum(arcs[inds, 3].^2 .* 0.5 .* (2 .* arcs[inds, 6] .- sin.(2 .* arcs[inds, 6])))
        total += _polygon_area(poly) + hole * circular
    end
    return total
end

"""Compute alpha-hull area; return `NaN` on invalid geometry like R's wrapper."""
function areaahull(h::AHull; timeout=5)
    try
        value = areaahulleval(h)
        return isfinite(value) && value >= 0 ? value : NaN
    catch
        return NaN
    end
end
