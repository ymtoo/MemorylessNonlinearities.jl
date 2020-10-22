module MemorylessNonlinearities

export Blanking, Cauchy, Clipping, HampelThreePart, TurkeyBiweight
export filt, minmaxrescale

include("utils.jl")

abstract type AbstractMemorylessNonlinearity end

struct Blanking{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct Cauchy{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
    isrescale::Bool
end
Cauchy(k) = Cauchy(k, false)

struct Clipping{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct HampelThreePart{T<:Real} <: AbstractMemorylessNonlinearity
    a::T
    b::T
    c::T
end

struct TurkeyBiweight{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

"""
Blanking nonlinearity.
"""
blanking(x::S, k::T) where {S<:Real,T<:Real} = abs(x) > k ? S(0) : x
function filt(f::Blanking, x::AbstractVector)
    blanking.(x, f.k)
end


"""
Cauchy nonlinearity.
"""
cauchy(x::S, k::T) where {S<:Real,T<:Real} = (2 * (x / k)) / (1 + (x / k)^2)
function filt(f::Cauchy, x::AbstractVector)
    y = cauchy.(x, f.k)
    f.isrescale ? minmaxrescale(y, -f.k, f.k) : y
end

"""
Clipping nonlinearity.
"""
clipping(x::S, k::T) where {S<:Real,T<:Real} = abs(x) > k ? sign(x) * k : x
function filt(f::Clipping, x::AbstractVector)
    clipping.(x, f.k)
end

"""
Hampel's Three Part Redescending nonlinearity.
"""
function hampelthreepart(x::S, a::T, b::T, c::T) where {S<:Real,T<:Real}
    absx = abs(x)
    if absx < a
        return x
    elseif absx >= a && absx < b
        return sign(x) * a 
    elseif absx >= b && absx < c
        return (absx - c) * (a / (b - c)) * sign(x)
    else
        return 0.
    end
end
function filt(f::HampelThreePart, x::AbstractVector)
    hampelthreepart.(x, f.a, f.b, f.c)
end

"""
Turkey's biweight nonlinearity.
"""
turkeybiweight(x::S, k::T) where {S<:Real,T<:Real} = abs(x) > k ? 0. : x * (1 - (x / k)^2)^2 
function filt(f::TurkeyBiweight, x::AbstractVector)
    turkeybiweight.(x, f.k)
end

end