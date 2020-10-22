using MemorylessNonlinearities
using Test



@testset "MemorylessNonlinearities.jl" begin

    x = [-3,-2,-1,0,1,2,3]
    @test filt(Blanking(1.0), x) == [0,0,-1,0,1,0,0]
    @test filt(Cauchy(1.0), x) == [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6]
    @test filt(Clipping(1.0), x) == [-1,-1,-1,0,1,1,1]
    @test filt(HampelThreePart(1.0,2.0,3.0), x) == [0.0,-1.0,-1.0,0.0,1.0,1.0,0.0]
    @test filt(TurkeyBiweight(3.0), x) ≈ [-0.0,-0.6173,-0.7901,0.0,0.7901,0.6173,0.0] atol=0.0001
    
end
