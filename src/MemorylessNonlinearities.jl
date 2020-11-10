module MemorylessNonlinearities

using AlphaStableDistributions, DelimitedFiles, Interpolations, QuadGK

export 
    Arctangent,
    Blanking, 
    CauchyNL, 
    Clipping, 
    HampelThreePart, 
    SαSNL, 
    SoftClipping, 
    TurkeyBiweight,

    filt, 
    minmaxrescale, 
    nlnames

include("utils.jl")

const DATADIR = joinpath(dirname(pathof(@__MODULE__)), "..", "data")

const _xsas = -200:0.1:200
const _αsas = 1.0:0.001:2.0
const _ysas = readdlm(joinpath(DATADIR, "sas.csv"), ',', Float64)

abstract type AbstractMemorylessNonlinearity end

struct Arctangent{T<:Real} <: AbstractMemorylessNonlinearity
    α::T
end

struct Blanking{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct CauchyNL{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct Clipping{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct HampelThreePart{T<:Real} <: AbstractMemorylessNonlinearity
    a::T
    b::T
    c::T
end

struct SαSNL{T<:Real} <: AbstractMemorylessNonlinearity 
    α::T
    scale::T
    location::T
    approx::Bool
end
SαSNL(α, scale, location) = SαSNL(α, scale, location, true)

struct SoftClipping{T<:Real} <: AbstractMemorylessNonlinearity 
    k::T
end

struct TurkeyBiweight{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

nlnames() = [:Arctangent,
             :Blanking, 
             :CauchyNL, 
             :Clipping, 
             :HampelThreePart, 
             :SαSNL, 
             :SoftClipping, 
             :TurkeyBiweight]

"""
Arctangent nonlinearity.
"""
arctangent(x::S, α::T) where {S<:Real,T<:Real} = (2 / π) * atan(α * x)
function filt(f::Arctangent, x::AbstractVector)
    arctangent.(x, f.α)
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
cauchynl(x::S, k::T) where {S<:Real,T<:Real} = (2 * (x / k)) / (1 + (x / k)^2)
function filt(f::CauchyNL, x::AbstractVector)
    cauchynl.(x, f.k)
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
Symmetric Alpha Stable nonlinearity.

`x` is standard symmetric alpha stable distributed.
"""
function sαsnl(x::S, α) where {S<:Real}
    num, _ = quadgk(t -> t * sin(t * x) * exp(-t^α), 0, Inf; atol=1e-12, order=3) 
    den, _ = quadgk(t -> cos(t * x) * exp(-t^α), 0, Inf; atol=1e-12, order=3)
    num / den
end

function filt(f::SαSNL, x::AbstractVector)
    (abs(f.α - 2.0) < 0.001) && return (x ./ 2) # Gaussian
    (abs(f.α - 1.0) < 0.001) && return (2 .* x) ./ (1 .+ x.^2) # Cauchy 
    xstd = (collect(x) .- f.location) ./ f.scale
    if f.approx
        maxval = maximum(_xsas)
        minval = minimum(_xsas)
        xstd[xstd .> maxval] .= (f.α+1) ./ xstd[xstd .> maxval]
        xstd[xstd .< minval] .= (f.α+1) ./ xstd[xstd .< minval]
        itp = interpolate((_αsas, _xsas), _ysas, Gridded(Linear()))
        xstd[(xstd .<= maxval) .& (xstd .>= minval)] .= itp(f.α, xstd[(xstd .<= maxval) .& (xstd .>= minval)])
        xstd
    else
        sαsnl.(xstd, f.α)
    end
end

"""
Soft clipping nonlinearity.
"""
function softclipping(x::S, k::T) where {S<:Real,T<:Real} 
    abs(x) > k ? sign(x) * k : k * (3*(x / k) / 2) * (1- ((x/ k)^2) / 3)
end
function filt(f::SoftClipping, x::AbstractVector)
    softclipping.(x, f.k)
end

"""
Turkey's biweight nonlinearity.
"""
turkeybiweight(x::S, k::T) where {S<:Real,T<:Real} = abs(x) > k ? 0. : x * (1 - (x / k)^2)^2 
function filt(f::TurkeyBiweight, x::AbstractVector)
    turkeybiweight.(x, f.k)
end

end