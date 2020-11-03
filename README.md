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

### Blanking 
![window](images/blanking-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, Blanking, ks)
```
![window](images/blanking.png)

### Cauchy
![window](images/cauchy-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, Cauchy, ks)
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

### SαS (approximated by 2D lookup table)
![window](images/sas-eqn.png)
```julia
αs = 1:0.2:2
plotnonlinearity(x, SαS, αs)
```
![window](images/sas.png)

### TurkeyBiweight
![window](images/turkeybiweight-eqn.png)
```julia
ks = 1:1:5
plotnonlinearity(x, TurkeyBiweight, ks)
```
![window](images/turkeybiweight.png)

## Performance
Chirp signals with Symmetric α-Stable noise parameterized by α=1.5, scale=1.0, location=0.0 are simulated. 

The following nonlinear functions are applied to filter the noise.
| Nonlinear       | Parameter                    |
| --------------- | ---------------------------- |
| Blanking        | k=3.0                        |
| Cauchy          | k=1.0                        |
| Clipping        | k=1.0                        |
| HampelThreePart | a=1.0,b=2.0,c=3.0            |
| SαS             | α=1.5,scale=1.0,location=0.0 |
| TurkeyBiweight  | k=3.0                        |

Root Mean Squared Errors (RMSEs) of filtered signals with respect to nonlinear functions and Generalizad Signal-to-Noise Ratios (GSNRs) are as follows. 
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