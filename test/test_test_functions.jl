let
    @test isapprox(ackley([0.0,0.0]), 0.0, atol=1e-14)
    @test isapprox(booth([1.0,3.0]),  0.0, atol=1e-14)
    @show branin(Float64[     -π,12.275])
    @show branin(Float64[      π, 2.275])
    @show branin(Float64[9.42478, 2.475])
    @test isapprox(branin(Float64[     -π,12.275]), 0.397887, atol=1e-6)
    @test isapprox(branin(Float64[      π, 2.275]), 0.397887, atol=1e-6)
    @test isapprox(branin(Float64[9.42478, 2.475]), 0.397887, atol=1e-6)
    @test isapprox(flower([0,0]), 0.0, atol=1e-14)
    @test isapprox(michalewicz([2.20319,1.57049]), -1.8013, atol=1e-4)
    @test isapprox(rosenbrock([1,1]), 0.0, atol=1e-14)
    @test wheeler([1.0,1.0])
    @test circle([0,1]) ≈ [0.0,1.0]
    @test circle([π,1]) ≈ [2.0,1.0]
end