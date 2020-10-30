module MemorylessNonlinearities

using AlphaStableDistributions, DelimitedFiles, Interpolations, QuadGK

export Blanking, Cauchy, Clipping, HampelThreePart, SαS, TurkeyBiweight
export filt, minmaxrescale

include("utils.jl")

const _xsas = -200:0.1:200
const _αsas = 1.0:0.001:2.0
const _ysas = readdlm(joinpath(@__DIR__, "..", "data/sas.csv"), ',', Float64)

abstract type AbstractMemorylessNonlinearity end

struct Blanking{T<:Real} <: AbstractMemorylessNonlinearity
    k::T
end

struct Cauchy{T<:Real} <: AbstractMemorylessNonlinearity
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

struct SαS{T<:Real} <: AbstractMemorylessNonlinearity 
    α::T
    scale::T
    location::T
    approx::Bool
end
SαS(α, scale, location) = SαS(α, scale, location, false)
SαS(α) = SαS(α, 1.0, 0.0, false)

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
    cauchy.(x, f.k)
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
function sαs(x::S, α) where {S<:Real}
    num, _ = quadgk(t -> t * sin(t * x) * exp(-t^α), 0, Inf; atol=1e-12, order=3) 
    den, _ = quadgk(t -> cos(t * x) * exp(-t^α), 0, Inf; atol=1e-12, order=3)
    num / den
end
K(α) = 13.0859 * α^4 - 68.4388 * α^3 + 134.7758 * α^2 - 115.9855 * α + 37.6752
function sαs(x::S, Kα, k) where {S<:Real} 
    τ = √(Kα/k)
    (abs(x) > τ) ? Kα / x : k*x
end

function filt(f::SαS, x::AbstractVector)
    (abs(f.α - 2.0) < 0.001) && return (x ./ 2) # Gaussian
    (abs(f.α - 1.0) < 0.001) && return (2 .* x) ./ (1 .+ x.^2) # Cauchy 
    xstd = (collect(x) .- f.location) ./ f.scale
    if f.approx
        # k = gamma(3/f.α) / gamma(1/f.α)
        # sαs.(xstd, Kα, k)
        maxval = maximum(_xsas)
        minval = minimum(_xsas)
        Kα = K(f.α)
        xstd[xstd .> maxval] .= (f.α+1) ./ xstd[xstd .> maxval]#Kα ./ xstd[xstd .> maxval]
        xstd[xstd .< minval] .= (f.α+1) ./ xstd[xstd .< minval]#Kα ./ xstd[xstd .< minval]
        itp = interpolate((_αsas, _xsas), _ysas, Gridded(Linear()))
        xstd[(xstd .<= maxval) .& (xstd .>= minval)] .= itp(f.α, xstd[(xstd .<= maxval) .& (xstd .>= minval)])
        xstd
    else
        sαs.(xstd, f.α)
    end
end

"""
Turkey's biweight nonlinearity.
"""
turkeybiweight(x::S, k::T) where {S<:Real,T<:Real} = abs(x) > k ? 0. : x * (1 - (x / k)^2)^2 
function filt(f::TurkeyBiweight, x::AbstractVector)
    turkeybiweight.(x, f.k)
end

end