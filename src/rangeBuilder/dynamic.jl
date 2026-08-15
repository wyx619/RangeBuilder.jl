const _range_world_cache = Dict{Int,Any}()
const _range_land_path = normpath(joinpath(@__DIR__, "..", "geodata", "ne_50m_land.jld2"))
# Match the R reference buffer default: 30 segments per quadrant.
# `getDynamicAlphaHull` calls its buffer operation without overriding
# nQuadSegs, so 30 segments per quadrant is the faithful equivalent.
# This affects only the polygonal approximation of a buffer.
const _R_BUFFER_QUADRANT_SEGMENTS = 30

mutable struct _RangeProjTransforms
    context::Ptr{Proj.PJ_CONTEXT}
    forward::Proj.Transformation
    inverse::Proj.Transformation
    geos::LibGEOS.GEOSContext
end

function _range_native_context()
    context = Proj.proj_context_create()
    forward = Proj.Transformation(
        "EPSG:4326", "+proj=eqearth";
        always_xy=true,
        ctx=context,
    )
    inverse = Proj.Transformation(
        "+proj=eqearth", "EPSG:4326";
        always_xy=true,
        ctx=context,
    )
    holder = _RangeProjTransforms(context, forward, inverse, LibGEOS.GEOSContext())
    finalizer(holder) do value
        finalize(value.forward)
        finalize(value.inverse)
        finalize(value.geos)
        value.context != C_NULL && Proj.proj_context_destroy(value.context)
        value.context = C_NULL
    end
    return holder
end

function _unique_finite_coordinates(points::AbstractMatrix)
    rows = Tuple{Float64,Float64}[]
    seen = Set{Tuple{Float64,Float64}}()
    for row in eachrow(points)
        all(isfinite, row) || continue
        point = (Float64(row[1]), Float64(row[2]))
        point in seen && continue
        push!(seen, point)
        push!(rows, point)
    end
    return isempty(rows) ? zeros(Float64, 0, 2) : reduce(vcat, ([p[1] p[2]] for p in rows))
end

function _range_shuffle_collinear(points::AbstractMatrix; rng=Random.MersenneTwister(0))
    # R rangeBuilder reshuffles all rows while the first three rows share an
    # exact x or y coordinate, so the Delaunay step is not collinear. R has no
    # loop bound; cap the attempts here to guarantee termination.
    size(points, 1) >= 3 || return points
    for _ in 1:100
        first_three = points[1:3, :]
        if first_three[1, 1] == first_three[2, 1] == first_three[3, 1] ||
           first_three[1, 2] == first_three[2, 2] == first_three[3, 2]
            points = points[Random.randperm(rng, size(points, 1)), :]
        else
            return points
        end
    end
    return points
end

function _closest_point_pair(points::AbstractMatrix)
    best_i, best_j, best_d2 = 1, 2, Inf
    for i in 1:size(points, 1), j in (i + 1):size(points, 1)
        dx = points[i, 1] - points[j, 1]
        dy = points[i, 2] - points[j, 2]
        d2 = dx * dx + dy * dy
        if d2 < best_d2
            best_i, best_j, best_d2 = i, j, d2
        end
    end
    return best_i, best_j
end

function _range_geodesic_closest_point_pair(points::AbstractMatrix)
    count = size(points, 1)
    unit_x = Vector{Float64}(undef, count)
    unit_y = Vector{Float64}(undef, count)
    unit_z = Vector{Float64}(undef, count)
    radians = pi / 180
    @inbounds for index in 1:count
        longitude = points[index, 1] * radians
        latitude = points[index, 2] * radians
        cos_latitude = cos(latitude)
        unit_x[index] = cos_latitude * cos(longitude)
        unit_y[index] = cos_latitude * sin(longitude)
        unit_z[index] = sin(latitude)
    end

    best_i, best_j, largest_dot = 1, 2, -2.0
    @inbounds for i in 1:(count - 1)
        xi, yi, zi = unit_x[i], unit_y[i], unit_z[i]
        for j in (i + 1):count
            dot = xi * unit_x[j] + yi * unit_y[j] + zi * unit_z[j]
            if dot > largest_dot
                best_i, best_j, largest_dot = i, j, dot
            end
        end
    end
    return best_i, best_j
end

function _range_drop_duplicate_points(points::AbstractMatrix; backend::Symbol=:delaunay)
    # Mirrors R rangeBuilder::getDynamicAlphaHull: while the Delaunay step
    # fails with a "duplicate data points" error, drop the closest point pair
    # and retry. Other failures (e.g. collinear input) are not retried here;
    # the caller falls back through its alpha loop exactly as R does.
    current = points
    while size(current, 1) >= 3
        if backend === :shull && !isnothing(_shull_float32_duplicate_pair(current))
            # R's rangeBuilder catches interp's "duplicate points" error,
            # then removes the first member of the geographically closest pair.
            # Detecting the Float32 collapse before SHull mutates its state
            # makes that retry deterministic and keeps the core pure.
            i, _ = _range_geodesic_closest_point_pair(current)
            current = current[[k for k in axes(current, 1) if k != i], :]
            continue
        end
        result = try
            delvor(current; backend)
        catch e
            e isa ArgumentError && occursin("duplicate data points", sprint(showerror, e)) ?
                nothing : return nothing
        end
        result !== nothing && return result
        i, j = _closest_point_pair(current)
        keep = [k for k in axes(current, 1) if k != i]
        current = current[keep, :]
    end
    return nothing
end

function _range_polygon_is_valid(poly, native=nothing)
    # Topological validity via the same GEOS backend used by the R reference.
    # Every component must be valid, matching R's `all(st_is_valid(hull))`.
    isnothing(poly) && return false
    components = _polygon_components(poly)
    isempty(components) && return false
    transforms = isnothing(native) ? _range_native_context() : native
    geos = transforms.geos
    for component in components
        trait = GeoInterface.geomtrait(component)
        geom = GeoInterface.convert(
            LibGEOS.geointerface_geomtype(trait), trait, component;
            context=geos,
        )
        LibGEOS.isValid(geom, geos) || return false
    end
    return true
end

function _range_projected_polygon_is_valid(poly, native=nothing)
    # In the dynamic-search loop, rangeBuilder transforms each newly-built
    # candidate to Equal Earth before calling st_is_valid() and st_buffer().
    # A valid buffer can cross the antimeridian after being transformed back
    # to EPSG:4326; GEOS then sees that geographic representation as a planar
    # self-intersection even though the projected geometry is valid.  Test the
    # candidate in the same working CRS as the R implementation.
    isnothing(poly) && return false
    components = _polygon_components(poly)
    isempty(components) && return false
    transforms = isnothing(native) ? _range_native_context() : native
    geos = transforms.geos
    for component in components
        projected = try
            GeometryOps.reproject(component, transforms.forward; target_crs="+proj=eqearth")
        catch
            return false
        end
        trait = GeoInterface.geomtrait(projected)
        projected_geos = try
            GeoInterface.convert(
                LibGEOS.geointerface_geomtype(trait), trait, projected;
                context=geos,
            )
        catch
            return false
        end
        try
            LibGEOS.isValid(projected_geos, geos) || return false
        catch
            return false
        end
    end
    return true
end

function _polygon_count(poly)
    isnothing(poly) && return 0
    trait = GeoInterface.geomtrait(poly)
    trait isa GeoInterface.MultiPolygonTrait && return GeoInterface.ngeom(trait, poly)
    trait isa GeoInterface.PolygonTrait && return 1
    trait isa GeoInterface.GeometryCollectionTrait && return GeoInterface.ngeom(trait, poly)
    return 0
end

@inline function _point_on_segment(px, py, ax, ay, bx, by)
    cross = (px - ax) * (by - ay) - (py - ay) * (bx - ax)
    abs(cross) <= 16eps(Float64) * max(1.0, abs(px), abs(py), abs(ax), abs(ay), abs(bx), abs(by)) || return false
    return min(ax, bx) <= px <= max(ax, bx) && min(ay, by) <= py <= max(ay, by)
end

function _ring_contains(px, py, ring)
    inside = false
    previous = last(ring)
    for current in ring
        ax, ay = Float64(previous[1]), Float64(previous[2])
        bx, by = Float64(current[1]), Float64(current[2])
        _point_on_segment(px, py, ax, ay, bx, by) && return true
        if (ay > py) != (by > py)
            x_intersection = ax + (py - ay) * (bx - ax) / (by - ay)
            x_intersection > px && (inside = !inside)
        end
        previous = current
    end
    return inside
end

function _polygon_contains(px, py, poly)
    exterior = GeoInterface.getexterior(poly)
    _ring_contains(px, py, GeoInterface.getpoint(exterior)) || return false
    for hole in GeoInterface.gethole(poly)
        _ring_contains(px, py, GeoInterface.getpoint(hole)) && return false
    end
    return true
end

function _geometry_contains(px, py, poly)
    trait = GeoInterface.geomtrait(poly)
    if trait isa GeoInterface.PolygonTrait
        return _polygon_contains(px, py, poly)
    elseif trait isa GeoInterface.MultiPolygonTrait || trait isa GeoInterface.GeometryCollectionTrait
        return any(_geometry_contains(px, py, GeoInterface.getgeom(trait, poly, i)) for i in 1:GeoInterface.ngeom(trait, poly))
    end
    return false
end

function _points_in_polygon(points::AbstractMatrix, poly, native=nothing)
    isnothing(poly) && return falses(size(points, 1))
    return [_geometry_contains(points[i, 1], points[i, 2], poly) for i in axes(points, 1)]
end

function _range_points_in_polygon(points::AbstractMatrix, poly, native=nothing)
    # Coverage predicate matching the R reference's planar intersection test
    # (sf_use_s2(FALSE)): a point intersects the polygon when it lies in its
    # interior or on any of its rings. The ray-casting fallback differs on
    # hole-boundary points, so the GEOS backend is authoritative.
    isnothing(poly) && return falses(size(points, 1))
    components = _polygon_components(poly)
    isempty(components) && return falses(size(points, 1))
    transforms = isnothing(native) ? _range_native_context() : native
    geos = transforms.geos
    results = falses(size(points, 1))
    for component in components
        ctrait = GeoInterface.geomtrait(component)
        cgeom = GeoInterface.convert(
            LibGEOS.geointerface_geomtype(ctrait), ctrait, component;
            context=geos,
        )
        for (idx, row) in enumerate(eachrow(points))
            point_geom = GeoInterface.convert(
                LibGEOS.Point, GeoInterface.PointTrait(),
                GeoInterface.Wrappers.Point((row[1], row[2])); context=geos,
            )
            results[idx] |= LibGEOS.intersects(point_geom, cgeom, geos)
        end
    end
    return results
end

function _range_point_intersection_counts(points::AbstractMatrix, poly, native=nothing)
    # `rangeBuilder` sums the per-point intersection list lengths.
    # Keep overlapping polygon features separate: a point intersecting two
    # components contributes two to coverage, rather than being collapsed to
    # a boolean membership result.
    isnothing(poly) && return zeros(Int, size(points, 1))
    components = _polygon_components(poly)
    isempty(components) && return zeros(Int, size(points, 1))
    transforms = isnothing(native) ? _range_native_context() : native
    geos = transforms.geos
    results = zeros(Int, size(points, 1))
    for component in components
        ctrait = GeoInterface.geomtrait(component)
        cgeom = GeoInterface.convert(
            LibGEOS.geointerface_geomtype(ctrait), ctrait, component;
            context=geos,
        )
        for (idx, row) in enumerate(eachrow(points))
            point_geom = GeoInterface.convert(
                LibGEOS.Point, GeoInterface.PointTrait(),
                GeoInterface.Wrappers.Point((row[1], row[2])); context=geos,
            )
            LibGEOS.intersects(point_geom, cgeom, geos) && (results[idx] += 1)
        end
    end
    return results
end

function _convex_hull_polygon(points::AbstractMatrix; crs=4326)
    size(points, 1) >= 3 || return nothing
    mesh_points = [Meshes.Point(points[i, 1], points[i, 2]) for i in axes(points, 1)]
    hull = Meshes.convexhull(mesh_points)
    ring = [(Float64(vertex.coords.x.val), Float64(vertex.coords.y.val))
            for vertex in Meshes.eachvertex(hull)]
    isempty(ring) && return nothing
    first(ring) == last(ring) || push!(ring, first(ring))
    return GeoInterface.Wrappers.Polygon([ring]; crs)
end

function _polygon_components(poly)
    trait = GeoInterface.geomtrait(poly)
    if trait isa GeoInterface.MultiPolygonTrait
        return [GeoInterface.getgeom(trait, poly, i) for i in 1:GeoInterface.ngeom(trait, poly)]
    end
    if trait isa GeoInterface.GeometryCollectionTrait
        return [GeoInterface.getgeom(trait, poly, i) for i in 1:GeoInterface.ngeom(trait, poly)]
    end
    trait isa GeoInterface.PolygonTrait && return [poly]
    return Any[]
end

function _buffer_range(poly, distance::Real; force::Bool=false, native=nothing)
    distance < 0 && throw(ArgumentError("buffer distance must be nonnegative"))
    components = _polygon_components(poly)
    isempty(components) && return poly
    # Do not short-circuit zero-width buffers. `rangeBuilder` always projects
    # candidates to Equal Earth and runs its buffer operation, including when
    # `buff = 0`. GEOS then normalizes the ring and the inverse projection can
    # move a point exactly on, or just outside, the boundary. That affects the
    # subsequent `st_intersects()` coverage check and may change alpha/MCH.
    # Each Julia task owns a PROJ context and reuses its transformations.
    # This avoids the Windows libproj global-context crash without serializing
    # candidate buffers.
    transforms = isnothing(native) ? _range_native_context() : native
    forward, inverse = transforms.forward, transforms.inverse
    geos = transforms.geos
    buffered = map(components) do component
        projected = GeometryOps.reproject(component, forward; target_crs="+proj=eqearth")
        trait = GeoInterface.geomtrait(projected)
        projected_geos = GeoInterface.convert(
            LibGEOS.geointerface_geomtype(trait), trait, projected;
            context=geos,
        )
        GeometryOps.reproject(
            LibGEOS.bufferWithStyle(
                projected_geos,
                Float64(distance);
                quadsegs=_R_BUFFER_QUADRANT_SEGMENTS,
                context=geos,
            ),
            inverse;
            target_crs="EPSG:4326"
        )
    end
    length(buffered) == 1 && return only(buffered)
    return GeoInterface.Wrappers.GeometryCollection(buffered; crs=4326)
end

function _naturalearth_land(scale::Integer)
    scale == 50 || throw(ArgumentError("only coastScale=50 is bundled"))
    haskey(_range_world_cache, scale) && return _range_world_cache[scale]
    isfile(_range_land_path) || error("bundled Natural Earth land data is missing: $_range_land_path")
    land_wkb = JLD2.load(_range_land_path, "land_wkb")
    world = LibGEOS.readgeom(Vector{Cuchar}(land_wkb))
    _range_world_cache[scale] = world
    return world
end

function _clip_range(poly, clipToCoast, scale::Integer)
    mode = clipToCoast isa Bool ? (clipToCoast ? :terrestrial : :no) : Symbol(clipToCoast)
    mode == :no && return poly
    mode in (:terrestrial, :aquatic) ||
        throw(ArgumentError("clipToCoast must be :no, :terrestrial, or :aquatic"))
    world = _naturalearth_land(scale)
    isnothing(world) && return poly
    return mode == :terrestrial ? GeometryOps.intersection(GeometryOps.GEOS(), poly, world) :
                                  GeometryOps.difference(GeometryOps.GEOS(), poly, world)
end

function _range_constraints(points::AbstractMatrix, poly, fraction::Real, partCount::Integer, native)
    coverage = sum(_range_point_intersection_counts(points, poly, native)) / size(points, 1)
    return _polygon_count(poly) <= partCount && coverage >= fraction
end

"""Generate a range polygon by increasing alpha until R-style constraints hold.

The return value is a named tuple `(hull, alpha)`, where `hull` is a
GeoInterface polygonal geometry and `alpha` is a string such as `"alpha7"`
or `"alphaMCH"`. When several original R-style polygon features are buffered
separately, `hull` is a geometry collection so their feature count is retained
for the dynamic `partCount` constraint.

The search mirrors `rangeBuilder::getDynamicAlphaHull` line by line: rows
whose first three points share an exact x or y are reshuffled, exact
duplicate points are dropped while the Delaunay step reports them, and alpha
advances until a topologically valid hull is produced. That hull is measured
for coverage before buffering; later candidates are buffered only when they
are valid, so coverage stays stale for invalid candidates exactly as in R.
"""
function getDynamicAlphaHull(x; fraction::Real=0.95, partCount::Integer=3,
                             buff::Real=10000, initialAlpha::Real=3,
                             coordHeaders=(1, 2), clipToCoast=:terrestrial,
                             alphaIncrement::Real=1, verbose::Bool=false,
                             alphaCap::Real=400, coastScale::Integer=50,
                             backend::Symbol=:delaunay)
    0 < fraction <= 1 || throw(ArgumentError("fraction must be in (0, 1]"))
    partCount >= 1 || throw(ArgumentError("partCount must be positive"))
    buff >= 0 || throw(ArgumentError("buff must be nonnegative"))
    alphaIncrement > 0 || throw(ArgumentError("alphaIncrement must be positive"))
    alphaCap >= initialAlpha || throw(ArgumentError("alphaCap must be at least initialAlpha"))

    points = _unique_finite_coordinates(_range_coordinates(x; coordHeaders))
    size(points, 1) >= 3 ||
        throw(ArgumentError("at least three unique finite coordinates are required"))

    # R reshuffles the rows while the first three share an exact x or y, then
    # drops duplicate points until the Delaunay step succeeds.
    points = _range_shuffle_collinear(points)
    delaunay = _range_drop_duplicate_points(points; backend)

    alpha = Float64(initialAlpha)
    hull = nothing
    accepted_alpha = nothing
    buffered = false
    cap = Float64(alphaCap) + eps(Float64(alphaCap))
    native = _range_native_context()

    alpha_polygon(alpha) = isnothing(delaunay) ? nothing : try
        ah2polygon(ahull(delaunay; alpha); crs=4326)
    catch
        nothing
    end

    # R L52-63/L90-93: advance alpha until a hull exists and is topologically
    # valid. The validity gate runs on the un-buffered candidate, matching
    # R's `ah2sf` + `st_is_valid` check before the first coverage measurement.
    hull = alpha_polygon(alpha)
    while alpha <= cap && (hull === nothing || !_range_polygon_is_valid(hull, native))
        verbose && println("\talpha: ", alpha)
        alpha += alphaIncrement
        alpha > cap && break
        hull = alpha_polygon(alpha)
    end

    if hull === nothing
        # R problem=TRUE: the alpha cap was exhausted, fall back to the
        # buffered convex hull of all points.
        hull = _buffer_range(_convex_hull_polygon(points; crs=4326), buff; force=true, native)
        accepted_alpha = "MCH"
        buffered = true
    else
        # R L97: initial coverage is measured on the un-buffered hull.
        pointWithin = _range_point_intersection_counts(points, hull, native)
        hull_valid = _range_polygon_is_valid(hull, native)

        # R L101: main loop while the feature count, coverage fraction, or
        # validity check fails. `hull_valid` records validity in the CRS where
        # the current hull was produced: geographic before the first buffer,
        # Equal Earth after a successful candidate buffer.
        while _polygon_count(hull) > partCount ||
              sum(pointWithin) / size(points, 1) < fraction ||
              !hull_valid
            alpha += alphaIncrement
            verbose && println("\talpha: ", alpha)
            candidate = alpha_polygon(alpha)

            # R L108-113: keep advancing while the candidate fails to build.
            while alpha <= cap && candidate === nothing
                alpha += alphaIncrement
                verbose && println("\talpha: ", alpha)
                candidate = alpha_polygon(alpha)
            end

            if candidate !== nothing
                # R L115-120: validate and buffer candidates in Equal Earth.
                # An inverse-projected buffer may cross the antimeridian, so
                # its EPSG:4326 coordinate representation is not used for the
                # next validity gate.
                if _range_projected_polygon_is_valid(candidate, native)
                    hull = _buffer_range(candidate, buff; force=true, native)
                    buffered = true
                    pointWithin = _range_point_intersection_counts(points, hull, native)
                    hull_valid = true
                else
                    hull = candidate
                    hull_valid = false
                end
            end

            # R L126-132: alpha past the cap selects the convex hull fallback.
            if alpha > cap
                hull = _buffer_range(_convex_hull_polygon(points; crs=4326), buff; force=true, native)
                accepted_alpha = "MCH"
                buffered = true
                hull_valid = true
                break
            end
        end

        accepted_alpha = something(accepted_alpha, alpha)

        # R L148-151: buffer the winning hull when it was never buffered.
        buffered || (hull = _buffer_range(hull, buff; force=true, native))
    end

    hull === nothing && throw(ArgumentError("could not construct a polygon from the coordinates"))
    hull = _clip_range(hull, clipToCoast, coastScale)
    alpha_label = accepted_alpha isa String ? accepted_alpha :
                  (isinteger(accepted_alpha) ? string(Int(round(accepted_alpha))) :
                   string(round(accepted_alpha; sigdigits=12)))
    return (hull=hull, alpha="alpha$(alpha_label)")
end
