include("../src/RootFinding.jl")

@static if VERSION < v"0.07-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "brent" begin
    f1(x) = x*x - 1
    f2(x) = x*x*x - x + 2
    f3(x) = sin(x)*exp(x*x-3.2) + 1/cos(x)
    @test f1(1) == 0
    @test f1(0) == -1
    res2 = -1.521379706804567569604081
    res3 = -3.14285155624251679480042160
    @test f2(res2) ≈ 0.0
    @test f2(0) != 0
    @test isapprox(RootFinding.brent(f1, -20.0, 0.0, 0.1, 0.1), -1, atol=0.1)
    @test isapprox(RootFinding.brent(f1, -2.0, 0.0, 0.1, 0.1), -1, atol=0.1)
    @test isapprox(RootFinding.brent(f1, -2.0, 0.0, 0.00001, 0.00001), -1, atol=0.00001)
    @test isapprox(RootFinding.brent(f1, -2.0, 0.0, 0.00001, 0.00001), -1, atol=0.00001)
    @test isapprox(RootFinding.brent(f1, -2.0, 0.0), -1)
    @test isapprox(RootFinding.rf(f1, -20.0, 0.0, 0.1, 0.1), -1, atol=0.1)
    @test isapprox(RootFinding.rf(f1, -2.0, 0.0, 0.1, 0.1), -1, atol=0.1)
    @test isapprox(RootFinding.rf(f1, -2.0, 0.0, 0.00001, 0.00001), -1, atol=0.00001)
    @test isapprox(RootFinding.rf(f1, -2.0, 0.0, 0.00001, 0.00001), -1, atol=0.00001)
    @test isapprox(RootFinding.rf(f1, -2.0, 0.0), -1)
    @test isapprox(RootFinding.brent(f2, -20.0, 0.0, 0.1, 0.1), res2, atol=0.1)
    @test isapprox(RootFinding.brent(f2, -2.0, 0.0, 0.1, 0.1), res2, atol=0.1)
    @test isapprox(RootFinding.brent(f2, -2.0, 0.0, 0.00001, 0.00001), res2, atol=0.00001)
    @test isapprox(RootFinding.brent(f2, -2.0, 0.0, 0.00001, 0.00001), res2, atol=0.00001)
    @test isapprox(RootFinding.brent(f2, -2.0, 0.0), res2)
    @test isapprox(RootFinding.rf(f2, -20.0, 0.0, 0.1, 0.1), res2, atol=0.1)
    @test isapprox(RootFinding.rf(f2, -2.0, 0.0, 0.1, 0.1), res2, atol=0.1)
    @test isapprox(RootFinding.rf(f2, -2.0, 0.0, 0.00001, 0.00001), res2, atol=0.00001)
    @test isapprox(RootFinding.rf(f2, -2.0, 0.0, 0.00001, 0.00001), res2, atol=0.00001)
    @test isapprox(RootFinding.rf(f2, -2.0, 0.0), res2)
    @test isapprox(RootFinding.brent(f3, -10.0, -2.1, 0.1, 0.1), res3, atol=0.1)
    @test isapprox(RootFinding.brent(f3, -10.0, -2.1, 0.1, 0.1), res3, atol=0.1)
    @test isapprox(RootFinding.brent(f3, -10.0, -2.1, 0.00001, 0.00001), res3, atol=0.00001)
    @test isapprox(RootFinding.brent(f3, -10.0, -2.1, 0.00001, 0.00001), res3, atol=0.00001)
    @test isapprox(RootFinding.brent(f3, -10.0, -2.1), res3)
    @test isapprox(RootFinding.rf(f3, -10.0, -2.1, 0.1, 0.1), -2.1, atol=0.1)
    @test isapprox(RootFinding.rf(f3, -10.0, -2.1, 0.1, 0.1), -2.1, atol=0.1)
    @test isapprox(RootFinding.rf(f3, -10.0, -2.1, 0.00001, 0.00001), -2.1, atol=0.00001)
    @test isapprox(RootFinding.rf(f3, -10.0, -2.1, 0.00001, 0.00001), -2.1, atol=0.00001)
    @test isapprox(RootFinding.rf(f3, -10.0, -2.1), -2.1)
end
