let
    x = [0.0]
    z = [1.0, 2.0]
    f = z -> sin(x[1]+z[1])*cos(x[1] + z[2])
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5])
    @test isapprox(ret[1], 0.0, atol=1e-10)
    @test isapprox(ret[2], 1.0, atol=1e-10)
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5], true)
    @test isapprox(ret[1], 0.0, atol=1e-10)
    @test isapprox(ret[2], 1.0, atol=1e-10)

    x = [-1.6]
    f = z -> sin(x[1]+z[1])*cos(x[1] + z[2])
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5])
    @test isapprox(ret[1], 0.0292, atol=1e-4)
    @test isapprox(ret[2], 0.4991, atol=1e-4)
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5], true)
    @test isapprox(ret[1], 0.0073, atol=1e-4)
    @test isapprox(ret[2], 0.5001, atol=1e-4)

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

    for (makepoly, prob, domain) in [
            (legendre, x->pdf(Uniform(-1,1),   x), (  -1,  1)),
            (laguerre, x->pdf(Exponential(1.0),x), (   0,Inf)),
            (hermite,  x->pdf(Normal(0.0,1.0), x), (-Inf,Inf)),
        ]
        bs = [Poly([1.0])]
        for i in 1 : 3
            push!(bs, orthogonal_recurrence(bs, prob, domain, 1e-16))
            b_pred = normalize(bs[end].a, 1)
            b_true = normalize(makepoly(i+1).a,1)
            @test min(norm(b_pred - b_true), norm(b_pred + b_true)) < 1e-10
        end
    end

    z = [1/3,2/3]
    bases = polynomial_chaos_bases([[legendre(1), legendre(2)], [laguerre(1), laguerre(2), laguerre(3)]])
    @test any([isapprox(b(z), legendre(1)(z[1])*laguerre(1)(z[2]), atol=1e-6) for b in bases])
    @test any([isapprox(b(z), legendre(1)(z[1])*laguerre(2)(z[2]), atol=1e-6) for b in bases])
    @test any([isapprox(b(z), legendre(1)(z[1])*laguerre(3)(z[2]), atol=1e-6) for b in bases])
    @test any([isapprox(b(z), legendre(2)(z[1])*laguerre(1)(z[2]), atol=1e-6) for b in bases])
    @test any([isapprox(b(z), legendre(2)(z[1])*laguerre(2)(z[2]), atol=1e-6) for b in bases])
    @test any([isapprox(b(z), legendre(2)(z[1])*laguerre(3)(z[2]), atol=1e-6) for b in bases])

    w = [1.0, 2.0]
    GP = GaussianProcess(x->1.0, (x,x′)->exp(-0.5*sum((xi-xi′)^2/wi for (xi,xi′,wi) in zip(x,x′,w))), [[1.0,1.0]], [3.0], 0.5)
    res = bayesian_monte_carlo(GP, w, [0.0,0.5], [1.0 -0.5; -0.5 2.0])
    @test isapprox(res[1], 1.2995, atol=1e-4)
    @test isapprox(res[2], 0.5457, atol=1e-4)
end