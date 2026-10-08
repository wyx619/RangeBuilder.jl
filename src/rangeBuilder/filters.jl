function _range_coordinates(x; coordHeaders=(1, 2))
    if x isa AbstractMatrix
        columns = collect(Int, coordHeaders)
        length(columns) == 2 || throw(ArgumentError("coordHeaders must select two columns"))
        all(1 .<= columns .<= size(x, 2)) || throw(ArgumentError("coordHeaders is out of bounds"))
        raw = x[:, columns]
        out = Matrix{Float64}(undef, size(raw, 1), 2)
        for i in axes(raw, 1), j in 1:2
            value = raw[i, j]
            out[i, j] = ismissing(value) ? NaN : Float64(value)
        end
        return out
    end

    headers = Symbol.(coordHeaders)
    length(headers) == 2 || throw(ArgumentError("coordHeaders must select two columns"))
    all(header -> hasproperty(x, header), headers) ||
        throw(ArgumentError("input does not provide the requested coordinate columns"))
    longitude = getproperty(x, headers[1])
    latitude = getproperty(x, headers[2])
    length(longitude) == length(latitude) || throw(ArgumentError("coordinate columns must have equal length"))
    out = Matrix{Float64}(undef, length(longitude), 2)
    for i in axes(out, 1)
        out[i, 1] = ismissing(longitude[i]) ? NaN : Float64(longitude[i])
        out[i, 2] = ismissing(latitude[i]) ? NaN : Float64(latitude[i])
    end
    return out
end

function _haversine_km(a::AbstractVector, b::AbstractVector)
    earth_radius_km = 6371.0088
    lat1, lon1 = deg2rad(a[2]), deg2rad(a[1])
    lat2, lon2 = deg2rad(b[2]), deg2rad(b[1])
    s_lat, s_lon = sin((lat2 - lat1) / 2), sin((lon2 - lon1) / 2)
    h = s_lat^2 + cos(lat1) * cos(lat2) * s_lon^2
    return 2earth_radius_km * asin(sqrt(clamp(h, 0.0, 1.0)))
end

function _is_lonlat_crs(crs)
    crs === 4326 && return true
    crs isa Integer && return crs == 4326
    value = lowercase(string(crs))
    return occursin("4326", value) || occursin("longlat", value) || occursin("wgs84", value)
end

"""Greedily remove points closer than `dist` kilometres.

With `returnIndex=true`, returns the one-based indices selected for removal;
when no points violate the threshold it returns `missing`, matching R
`filterByProximity`.  Matrix inputs are interpreted as longitude/latitude by
default, while projected coordinate matrices use metres when `crs` is changed.
"""
function filterByProximity(xy, dist::Real; returnIndex::Bool=false,
                           coordHeaders=(1, 2), crs=4326)
    dist >= 0 || throw(ArgumentError("dist must be nonnegative"))
    points = _range_coordinates(xy; coordHeaders)
    all(isfinite, points) || throw(ArgumentError("coordinates must be finite"))
    n = size(points, 1)
    n == 0 && return returnIndex ? missing : points

    pairs = Tuple{Int,Int}[]
    lonlat = _is_lonlat_crs(crs)
    for i in 1:n
        for j in 1:n
            i == j && continue
            separation = lonlat ? _haversine_km(view(points, i, :), view(points, j, :)) :
                                  hypot(points[i, 1] - points[j, 1], points[i, 2] - points[j, 2]) / 1000
            separation <= dist && push!(pairs, (i, j))
        end
    end
    isempty(pairs) && return returnIndex ? missing : points

    discard = Int[]
    while !isempty(pairs)
        point = pairs[1][1]
        push!(discard, point)
        filter!(pair -> pair[1] != point && pair[2] != point, pairs)
    end
    return returnIndex ? discard : points[setdiff(1:n, discard), :]
end
