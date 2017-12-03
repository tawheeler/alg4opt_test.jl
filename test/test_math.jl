let
    quadrule  = quadrule_legendre(3)
    i = findfirst(quadrule.xs, 0.0)
    @test isapprox(quadrule.ws[i], 8/9, atol=1e-8)
    i = findfirst(x->isapprox(x,  sqrt(3/5), atol=1e-6), quadrule.xs)
    @test isapprox(quadrule.ws[i], 5/9, atol=1e-8)
    i = findfirst(x->isapprox(x, -sqrt(3/5), atol=1e-6), quadrule.xs)
    @test isapprox(quadrule.ws[i], 5/9, atol=1e-8)

    f = x -> x^5- 2x^4 + 3x^3 + 5x^2 -x + 4
    @test isapprox(quadint(f, quadrule_legendre(3)), 10.533, atol=1e-3)
    @test isapprox(quadint(f, quadrule_legendre(3), -3, 5), 1820.8, atol=1e-3)
end