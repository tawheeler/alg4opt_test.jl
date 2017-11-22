let
    @test legendre(1) == Poly([1])
    @test legendre(2) == Poly([0,1])
    @test legendre(3) == Poly([-1/2,0,3/2])
    @test legendre(4) == Poly([0,-3/2,0,5/2])

    @test laguerre(1) == Poly([1])
    @test laguerre(2) == Poly([1,-1])
    @test laguerre(3) == Poly([1,-2,1/2])
    @test laguerre(4) == Poly([1,-3,3/2,-1/6])

    @test hermite(1) == Poly([1])
    @test hermite(2) == Poly([0,1])
    @test hermite(3) == Poly([-1,0,1])
    @test hermite(4) == Poly([0,-3,0,1])

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