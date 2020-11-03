using MemorylessNonlinearities
using Test

x = [-3,-2,-1,0,1,2,3]

@testset "MemorylessNonlinearities.jl" begin

    @info "Test MemorylessNonlinearities"

    @test filt(Blanking(1.0), x) == [0,0,-1,0,1,0,0]

    @test filt(Cauchy(1.0), x) == [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6]

    @test filt(Clipping(1.0), x) == [-1,-1,-1,0,1,1,1]

    @test filt(HampelThreePart(1.0,2.0,3.0), x) == [0.0,-1.0,-1.0,0.0,1.0,1.0,0.0]

    @test filt(SαS(2.0, 1.0, 0.0, false), x) ≈ x ./ 2 atol=0.0001
    @test filt(SαS(1.0, 1.0, 0.0, false), x) ≈ [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6] atol=0.0001
    @test filt(SαS(2.0, 1.0, 0.0, true), x) ≈ x ./ 2 atol=0.0001
    @test filt(SαS(1.0, 1.0, 0.0, true), x) ≈ [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6] atol=0.0001
    @test filt(SαS(2.0, 1.0, 0.0), x) ≈ x ./ 2 atol=0.0001
    @test filt(SαS(1.0, 1.0, 0.0), x) ≈ [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6] atol=0.0001
    @test filt(SαS(2.0), x) ≈ x ./ 2 atol=0.0001
    @test filt(SαS(1.0), x) ≈ [-0.6,-0.8,-1.0,0.0,1.0,0.8,0.6] atol=0.0001
    αs = [1.000001, 1.001000, 1.002000, 1.123457, 1.234567, 1.345678, 1.456789, 
          1.567899, 1.678999, 1.789999, 1.899999, 1.998000, 1.999000, 1.999999]
    xtmp = [-999.987654, -500.505050, -100.101010, -30.303030, -10.101010, -5.505050, -3.303030, -1.101010, -0.010101,
            0.0, 0.010101, 1.101010, 3.303030, 5.505050, 10.101010, 30.303030, 100.101010, 500.505050, 999.987654]
    for α in αs
        @test filt(SαS(α, 1.0, 0.0, true), xtmp) ≈ filt(SαS(α, 1.0, 0.0, false), xtmp) rtol=0.01
        @test filt(SαS(α, 2.0, 0.0), xtmp) ≈ filt(SαS(α, 2.0, 0.0, false), xtmp) rtol=0.01
        @test filt(SαS(α), xtmp) ≈ filt(SαS(α, 1.0, 0.0, false), xtmp) rtol=0.01
    end

    @test filt(TurkeyBiweight(3.0), x) ≈ [-0.0,-0.6173,-0.7901,0.0,0.7901,0.6173,0.0] atol=0.0001
    
end

@testset "utils.jl" begin

    @info "Test utils"

    x = randn(1000)
    xrescale = minmaxrescale(x, -0.1, 0.3)
    @test minimum(xrescale) ≈ -0.1 atol=1e-6
    @test maximum(xrescale) ≈ 0.3 atol=1e-6

end
