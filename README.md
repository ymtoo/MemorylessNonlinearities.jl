# MemorylessNonlinearities

[![Build Status](https://travis-ci.com/ymtoo/MemorylessNonlinearities.jl.svg?branch=master)](https://travis-ci.com/ymtoo/MemorylessNonlinearities.jl)
[![Coverage](https://codecov.io/gh/ymtoo/MemorylessNonlinearities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ymtoo/MemorylessNonlinearities.jl)

This package implements memoryless nonlinearity functions.

## Usage

```julia
using MemorylessNonlinearities, Plots

x = -10:0.1:10

function plotnonlinearity(x, f, params; size=(800, 600))
    p = plot(size=size, legend=:outertopright)
    for param in params
        label = join([String(name) * "=$(p)" for (name, p) in zip(
                     fieldnames(f), param)], ",")
        plot!(p, x, minmaxrescale(filt(f(param...), x), -1.0, 1.0); 
              linewidth=2, label=label)
    end
    p
end
```

### Arctangent
![window](images/arctangent-eqn.png)
```julia
αs = 1:1:10
plotnonlinearity(x, Arctangent, αs)
```
![window](images/arctangent.png)

### Blanking 
![window](images/blanking-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, Blanking, ks)
```
![window](images/blanking.png)

### CauchyNL
![window](images/cauchy-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, CauchyNL, ks)
```
![window](images/cauchy.png)

### Clipping
![window](images/clipping-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, Clipping, ks)
```
![window](images/clipping.png)

### HampelThreePart
![window](images/hampelthreepart-eqn.png)
```julia
abcs = ((1, 2, 3), (2, 3, 4), (3, 4, 5))
plotnonlinearity(x, HampelThreePart, abcs)
```
![window](images/hampelthreepart.png)

### SαSNL
Nonlinearity in locally optimal detectors based on IID SαS noise(approximated by 2D lookup table)
![window](images/sas-eqn.png)
```julia
αs = 1:0.2:2
plotnonlinearity(x, SαSNL, αs)
```
![window](images/sas.png)

### SoftClipping
![window](images/softclipping-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, SoftClipping, ks)
```
![window](images/softclipping.png)

### TurkeyBiweight
![window](images/turkeybiweight-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, TurkeyBiweight, ks)
```
![window](images/turkeybiweight.png)

## Performance
Chirp signals with Symmetric α-Stable noise parameterized by α=1.5, scale=1.0, location=0.0 were simulated. The following nonlinear functions were applied to the simulated data to filter the noise.
| Nonlinear       | Parameter                 |
| --------------- | ------------------------- |
| Arctangent      | α=1 
| Blanking        | k=3σ                      |
| CauchyNL        | k=3σ                      |
| Clipping        | k=σ                       |
| HampelThreePart | a=3σ,b=4σ,c=5σ            |
| SαSNL           | α=α',scale=c',location=δ' |
| SoftClipping    | k=σ                       |
| TurkeyBiweight  | k=3σ                      |

σ is median absolution deviation of the simulated data. a', c' and δ' are the estimated pararamters of IID Symmetric α-Stable distributions based on the simulated data. Root Mean Squared Errors (RMSEs) between the true chirp signals and filtered signals with respect to nonlinear functions and Generalizad Signal-to-Noise Ratios (GSNRs) are as follows. 
```julia
include("perf/simulate.jl")

E, gsnrs, gs = simulate()
plot(gsnrs, 
     dropdims(sum(E, dims=1) / size(E, 1), dims=1); 
     linewidth=2,
     label=reshape(string.(first.(gs)), 1, length(gs)),
     size=(800, 600), 
     legend=:outertopright, 
     xlabel="GSNR",
     ylabel="RMSE")
```
![window](images/rmse.png)