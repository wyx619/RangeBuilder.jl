const _range_world_cache = Dict{Int,Any}()
const _range_land_path = normpath(joinpath(@__DIR__, "..", "geodata", "ne_50m_land.jld2"))
# Match terra::buffer(), which approximates circular arcs with ten segments
# per quadrant. This affects only the polygonal approximation of a buffer.
const _R_BUFFER_QUADRANT_SEGMENTS = 10

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
    # A zero-width buffer is the identity operation. In particular, do not
    # project it: R's `buff=0` path does not change candidate geometry.
    distance == 0 && return poly
    distance < 0 && throw(ArgumentError("buffer distance must be nonnegative"))
    components = _polygon_components(poly)
    isempty(components) && return poly
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
    coverage = count(_points_in_polygon(points, poly, native)) / size(points, 1)
    return _polygon_count(poly) <= partCount && coverage >= fraction
end

"""Generate a range polygon by increasing alpha until R-style constraints hold.

The return value is a named tuple `(hull, alpha)`, where `hull` is a
GeoInterface polygonal geometry and `alpha` is a string such as `"alpha7"`
or `"alphaMCH"`. When several original R-style polygon features are buffered
separately, `hull` is a geometry collection so their feature count is retained
for the dynamic `partCount` constraint.
"""
function getDynamicAlphaHull(x; fraction::Real=0.95, partCount::Integer=3,
                             buff::Real=10000, initialAlpha::Real=3,
                             coordHeaders=(1, 2), clipToCoast=:terrestrial,
                             alphaIncrement::Real=1, verbose::Bool=false,
                             alphaCap::Real=400, coastScale::Integer=50)
    0 < fraction <= 1 || throw(ArgumentError("fraction must be in (0, 1]"))
    partCount >= 1 || throw(ArgumentError("partCount must be positive"))
    buff >= 0 || throw(ArgumentError("buff must be nonnegative"))
    alphaIncrement > 0 || throw(ArgumentError("alphaIncrement must be positive"))
    alphaCap >= initialAlpha || throw(ArgumentError("alphaCap must be at least initialAlpha"))

    points = _unique_finite_coordinates(_range_coordinates(x; coordHeaders))
    size(points, 1) >= 3 ||
        throw(ArgumentError("at least three unique finite coordinates are required"))

    alpha = Float64(initialAlpha)
    hull = nothing
    accepted_alpha = nothing
    buffered = false
    cap = Float64(alphaCap) + eps(Float64(alphaCap))
    native = _range_native_context()
    delaunay = try
        delvor(points)
    catch
        nothing
    end

    alpha_polygon(alpha) = isnothing(delaunay) ? nothing : try
        ah2polygon(ahull(delaunay; alpha); crs=4326)
    catch
        nothing
    end

    # R evaluates its first valid alpha hull before buffering it. Once that
    # candidate fails the constraints, every later candidate is buffered in
    # Equal Earth coordinates before its coverage is measured.
    candidate = nothing
    while alpha <= cap && candidate === nothing
        verbose && println("\talpha: ", alpha)
        candidate = alpha_polygon(alpha)
        candidate === nothing && (alpha += alphaIncrement)
    end

    if candidate !== nothing && _range_constraints(points, candidate, fraction, partCount, native)
        hull = candidate
        accepted_alpha = alpha
    end

    while hull === nothing && alpha <= cap
        alpha += alphaIncrement
        alpha > cap && break
        verbose && println("\talpha: ", alpha)
        candidate = alpha_polygon(alpha)
        candidate === nothing && continue
        candidate = _buffer_range(candidate, buff; force=true, native)
        if _range_constraints(points, candidate, fraction, partCount, native)
            hull = candidate
            accepted_alpha = alpha
            buffered = true
        end
    end

    if hull === nothing
        hull = _convex_hull_polygon(points; crs=4326)
        accepted_alpha = "MCH"
        buffered = false
    end
    hull === nothing && throw(ArgumentError("could not construct a polygon from the coordinates"))
    buffered || (hull = _buffer_range(hull, buff; force=true, native))
    hull = _clip_range(hull, clipToCoast, coastScale)
    alpha_label = accepted_alpha isa String ? accepted_alpha :
                  (isinteger(accepted_alpha) ? string(Int(round(accepted_alpha))) :
                   string(round(accepted_alpha; sigdigits=12)))
    return (hull=hull, alpha="alpha$(alpha_label)")
end
