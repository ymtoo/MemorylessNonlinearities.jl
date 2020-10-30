# MemorylessNonlinearities

[![Build Status](https://travis-ci.com/ymtoo/MemorylessNonlinearities.jl.svg?branch=master)](https://travis-ci.com/ymtoo/MemorylessNonlinearities.jl)
[![Coverage](https://codecov.io/gh/ymtoo/MemorylessNonlinearities.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ymtoo/MemorylessNonlinearities.jl)

This package implements memoryless nonlinearity functions.

## Usage
```julia
using MemorylessNonlinearities, Plots

x = -10:0.1:10
gs = [(Blanking, 3.0),
      (Cauchy, 3.0),
      (Clipping, 3.0),
      (HampelThreePart, (3.0, 4.0, 5.0)),
      (SαS, 1.5),
      (TurkeyBiweight, 5.0)]
p = plot(size=(800, 600), legend=:outertopright)
for (g, params) in gs
    plot!(p, x, minmaxrescale(filt(g(params...), x), -1.0, 1.0); 
          linewidth=2, label=string(g))
end
p
```
![window](nonlinearities.png)
```julia
x = -10:0.1:10
αs = 1.0:0.1:2.0
p = plot(size=(800, 600), legend=:outertopright)
for α in αs
    plot!(p, x, minmaxrescale(filt(SαS(α), x), -1.0, 1.0); 
          linewidth=2, label="α=$(α)")
end
p
```
![window](sas.png)