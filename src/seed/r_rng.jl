module RSeed

export AbstractRUniformRNG, RMersenneTwister, RLecuyerCMRG, r_rng, r_unif, r_runif, r_jitter

const _N = 624
const _M = 397
const _MATRIX_A = UInt32(0x9908b0df)
const _UPPER_MASK = UInt32(0x80000000)
const _LOWER_MASK = UInt32(0x7fffffff)
const _TEMPERING_MASK_B = UInt32(0x9d2c5680)
const _TEMPERING_MASK_C = UInt32(0xefc60000)
const _I2_32M1 = 2.328306437080797e-10
const _MT_SCALE = 2.3283064365386963e-10

abstract type AbstractRUniformRNG end

"""State of R 4.6's default `Mersenne-Twister` uniform generator."""
mutable struct RMersenneTwister <: AbstractRUniformRNG
    state::Vector{UInt32}
    index::Int
end

"""State of R 4.6's `L'Ecuyer-CMRG` uniform generator."""
mutable struct RLecuyerCMRG <: AbstractRUniformRNG
    state::Vector{UInt32}
end

@inline _lcg(value::UInt32) = UInt32(69069) * value + UInt32(1)

function RMersenneTwister(seed::Integer)
    value = UInt32(mod(Int128(seed), Int128(1) << 32))
    for _ in 1:50
        value = _lcg(value)
    end
    # R stores i_seed[0] first and then the 624 MT state words. FixupSeeds()
    # replaces i_seed[0] with 624 before the first generation cycle.
    value = _lcg(value)
    state = Vector{UInt32}(undef, _N)
    @inbounds for index in 1:_N
        value = _lcg(value)
        state[index] = value
    end
    return RMersenneTwister(state, _N)
end

"""Construct an R-compatible generator from R's full `.Random.seed` vector.

Only the Mersenne-Twister variant is supported. The vector must be an
unmodified R state with its kind code, current state index, and 624 signed
32-bit state words.
"""
function RMersenneTwister(r_seed::AbstractVector{<:Integer})
    length(r_seed) == _N + 2 ||
        throw(ArgumentError("an R Mersenne-Twister .Random.seed must contain $(_N + 2) integers"))
    mod(r_seed[1], 100) == 3 ||
        throw(ArgumentError(".Random.seed does not use R's Mersenne-Twister generator"))
    0 <= r_seed[2] <= _N ||
        throw(ArgumentError("invalid R Mersenne-Twister state index"))

    state = Vector{UInt32}(undef, _N)
    @inbounds for index in eachindex(state)
        state[index] = reinterpret(UInt32, Int32(r_seed[index + 2]))
    end
    return RMersenneTwister(state, Int(r_seed[2]))
end

const _LECUYER_M1 = Int64(4_294_967_087)
const _LECUYER_M2 = Int64(4_294_944_443)
const _LECUYER_NORM = 2.328306549295727688e-10
const _LECUYER_A12 = Int64(1_403_580)
const _LECUYER_A13N = Int64(810_728)
const _LECUYER_A21 = Int64(527_612)
const _LECUYER_A23N = Int64(1_370_589)

function RLecuyerCMRG(seed::Integer)
    value = UInt32(mod(Int128(seed), Int128(1) << 32))
    for _ in 1:50
        value = _lcg(value)
    end
    state = Vector{UInt32}(undef, 6)
    @inbounds for index in eachindex(state)
        value = _lcg(value)
        while Int64(value) >= _LECUYER_M2
            value = _lcg(value)
        end
        state[index] = value
    end
    return RLecuyerCMRG(state)
end

"""Construct an `L'Ecuyer-CMRG` generator from R's full `.Random.seed` vector."""
function RLecuyerCMRG(r_seed::AbstractVector{<:Integer})
    length(r_seed) == 7 ||
        throw(ArgumentError("an R L'Ecuyer-CMRG .Random.seed must contain 7 integers"))
    mod(r_seed[1], 100) == 7 ||
        throw(ArgumentError(".Random.seed does not use R's L'Ecuyer-CMRG generator"))
    state = Vector{UInt32}(undef, 6)
    @inbounds for index in eachindex(state)
        state[index] = reinterpret(UInt32, Int32(r_seed[index + 1]))
    end
    all(Int64(value) < _LECUYER_M1 for value in state[1:3]) ||
        throw(ArgumentError("invalid first L'Ecuyer-CMRG state component"))
    all(Int64(value) < _LECUYER_M2 for value in state[4:6]) ||
        throw(ArgumentError("invalid second L'Ecuyer-CMRG state component"))
    return RLecuyerCMRG(state)
end

"""Construct the supported R uniform generator encoded by `.Random.seed`."""
function r_rng(r_seed::AbstractVector{<:Integer})
    isempty(r_seed) && throw(ArgumentError(".Random.seed must not be empty"))
    kind = mod(r_seed[1], 100)
    kind == 3 && return RMersenneTwister(r_seed)
    kind == 7 && return RLecuyerCMRG(r_seed)
    throw(ArgumentError("unsupported R uniform generator kind $kind"))
end

@inline function _fixup(value::Float64)
    value <= 0.0 && return 0.5 * _I2_32M1
    (1.0 - value) <= 0.0 && return 1.0 - 0.5 * _I2_32M1
    return value
end

function _twist!(rng::RMersenneTwister)
    state = rng.state
    @inbounds for index in 1:(_N - _M)
        y = (state[index] & _UPPER_MASK) | (state[index + 1] & _LOWER_MASK)
        state[index] = xor(state[index + _M], y >> 1, (y & UInt32(1)) == 0 ? UInt32(0) : _MATRIX_A)
    end
    @inbounds for index in (_N - _M + 1):(_N - 1)
        y = (state[index] & _UPPER_MASK) | (state[index + 1] & _LOWER_MASK)
        state[index] = xor(state[index + (_M - _N)], y >> 1, (y & UInt32(1)) == 0 ? UInt32(0) : _MATRIX_A)
    end
    @inbounds begin
        y = (state[_N] & _UPPER_MASK) | (state[1] & _LOWER_MASK)
        state[_N] = xor(state[_M], y >> 1, (y & UInt32(1)) == 0 ? UInt32(0) : _MATRIX_A)
    end
    rng.index = 0
    return nothing
end

function _mt_genrand!(rng::RMersenneTwister)
    rng.index >= _N && _twist!(rng)
    rng.index += 1
    y = @inbounds rng.state[rng.index]
    y = xor(y, y >> 11)
    y = xor(y, (y << 7) & _TEMPERING_MASK_B)
    y = xor(y, (y << 15) & _TEMPERING_MASK_C)
    y = xor(y, y >> 18)
    return Float64(y) * _MT_SCALE
end

function _lecuyer_unif!(rng::RLecuyerCMRG)
    state = rng.state
    @inbounds begin
        p1 = _LECUYER_A12 * Int64(state[2]) - _LECUYER_A13N * Int64(state[1])
        p1 -= div(p1, _LECUYER_M1) * _LECUYER_M1
        p1 < 0 && (p1 += _LECUYER_M1)
        state[1], state[2], state[3] = state[2], state[3], UInt32(p1)

        p2 = _LECUYER_A21 * Int64(state[6]) - _LECUYER_A23N * Int64(state[4])
        p2 -= div(p2, _LECUYER_M2) * _LECUYER_M2
        p2 < 0 && (p2 += _LECUYER_M2)
        state[4], state[5], state[6] = state[5], state[6], UInt32(p2)

        return Float64(p1 > p2 ? p1 - p2 : p1 - p2 + _LECUYER_M1) * _LECUYER_NORM
    end
end

"""Generate one R-compatible open-interval uniform variate."""
function r_unif(rng::RMersenneTwister)
    return _fixup(_mt_genrand!(rng))
end

function r_unif(rng::RLecuyerCMRG)
    return _lecuyer_unif!(rng)
end

"""Generate `n` values matching R's `runif(n, min, max)`."""
function r_runif(rng::AbstractRUniformRNG, n::Integer, min::Real=0.0, max::Real=1.0)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    lower, upper = Float64(min), Float64(max)
    isfinite(lower) && isfinite(upper) && upper >= lower ||
        throw(ArgumentError("invalid uniform bounds"))
    lower == upper && return fill(lower, n)
    result = Vector{Float64}(undef, n)
    @inbounds for index in eachindex(result)
        value = r_unif(rng)
        result[index] = lower + (upper - lower) * value
    end
    return result
end

"""Apply R's `jitter` transform with an already computed amount."""
function r_jitter(rng::AbstractRUniformRNG, values::AbstractVector{<:Real}, amount::Real)
    offsets = r_runif(rng, length(values), -1.0, 1.0)
    return Float64.(values) .+ Float64(amount) .* offsets
end

end
