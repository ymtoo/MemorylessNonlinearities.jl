using AlphaStableDistributions
using Distributions
using DSP
using MemorylessNonlinearities
using SignalAnalysis
using StatsBase

export simulate

const _α = 1.5
const _scale = 1.
const _location = 0.

function rand_chirp(fs, recordtlen, chirptlen, freqs::Tuple)
    N = trunc(Integer, fs * recordtlen)
    M = trunc(Integer, fs * chirptlen)
    x = zeros(N)
    startind = N-M == 0 ? 1 : rand(1:N-M)
    startfreq = rand(Uniform(freqs...))
    stopfreq = rand(Uniform(startfreq, freqs[end]))
    x[startind:startind+M-1] = real(samples(chirp(startfreq, stopfreq, chirptlen, fs)))
    x
end

function get_attenuation(signal::AbstractVector{T}, noise::AbstractVector{T}, gsnr) where {T<:Real}
    Cg = 1.7811
    rmssignal = sqrt(sum(abs2, signal) / length(signal))
    So = exp(sum(log.(eps(eltype(noise)) .+ abs.(noise))) / length(noise))
    A = sqrt(4 * Cg * (So^2) * 10^(gsnr/10))
    A / rmssignal
end

function simulate(gs=nothing;
                  nrealizations=50,
                  gsnrs = -20:1:20,
                  recordtlen=1, 
                  chirptlen=1, 
                  freqs=(100, 2000),
                  fs=9600)
    # (gs === nothing) && (gs = [(Blanking, 3 * _scale), 
    #                            (CauchyNL, 3 * _scale), 
    #                            (Clipping, 3 * _scale),
    #                            (HampelThreePart, 3 * _scale, 6 * _scale, 9 * _scale),
    #                            (SαSNL, _α, 3 * _scale, _location),
    #                            (TurkeyBiweight, 3 * _scale)])
    ngsnrs = length(gsnrs)
    ngs = gs === nothing ? length(nlnames()) : length(gs)
    E = zeros(nrealizations, ngsnrs, ngs)
    for i in 1:ngsnrs
        for j in 1:nrealizations
            s = rand_chirp(fs, recordtlen, chirptlen, freqs)
            v = rand(AlphaStable(α=_α, scale=_scale, location=_location), length(s))
            atten = get_attenuation(s, v, gsnrs[i])
            x = atten .* s .+ v
            σ = mad(x)
            d = fit(SymmetricAlphaStable, x)
            (gs === nothing) && (gs = [(Blanking, 3 * σ), 
                                 (CauchyNL, 3 * σ), 
                                 (Clipping, σ),
                                 (HampelThreePart, 3 * σ, 4 * σ, 5 * σ),
                                 (SαSNL, d.α, d.scale, d.location),
                                 (TurkeyBiweight, 3 * σ)])
            for (k, g) in enumerate(gs)
                ŝ = MemorylessNonlinearities.filt(g[1](g[2:end]...), x)
                a = (ŝ's) / (ŝ'ŝ)
                E[j, i, k] = rms(s - a .* ŝ)
            end
        end
    end
    E, gsnrs, gs
end