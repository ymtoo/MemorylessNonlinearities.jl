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

```julia
ks = 1:1:5
plotnonlinearity(x, Blanking, ks)
```
![window](images/blanking.png)

### Cauchy
```julia
ks = 1:1:5
plotnonlinearity(x, Cauchy, ks)
```
![window](images/cauchy.png)

### Clipping
```julia
ks = 1:1:5
plotnonlinearity(x, Clipping, ks)
```
![window](images/clipping.png)

### HampelThreePart
```julia
abcs = ((1, 2, 3), (2, 3, 4), (3, 4, 5))
plotnonlinearity(x, HampelThreePart, abcs)
```
![window](images/hampelthreepart.png)

### SαS (approximation)
```julia
αs = 1:0.2:2
plotnonlinearity(x, SαS, αs)
```
![window](images/sas.png)

### TurkeyBiweight
```julia
ks = 1:1:5
plotnonlinearity(x, TurkeyBiweight, ks)
```
![window](images/turkeybiweight.png)