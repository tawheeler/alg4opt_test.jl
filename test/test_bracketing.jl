let
    a, b = bracket_minimum(x->(x-1)^2)
    @test a < 1 < b

    a, b = fibonacci_search(x->x^2, 0, 1, 2)
    @test a == 0
    @test b ≈ 0.5
    a, b = fibonacci_search(x->x^2, 0, 3, 3)
    @test a == 0
    @test b ≈ 1.0
    a, b = fibonacci_search(x->x^2, 0, 5, 4)
    @test a == 0
    @test b ≈ 1.0

    a, b = golden_section_search(x->x^2, 0, φ, 2)
    @test a == 0
    @test b ≈ 1
    a, b = golden_section_search(x->x^2, 0, φ^2, 3)
    @test a == 0
    @test b ≈ 1

    a, b, c = quadratic_fit_search(x->x^2, -1, -0.25, 1, 4)
    @test a == -0.25
    @test b ≈ 0.0
    @test c == 1.0

    P, intervals = shubert_piyavskii(x->sin(x) + sin(10/3*x), -2.7, 7.5, 6.5, 0.01)
    @test isapprox(P.x,  5.1452, atol=1e-3)
    @test isapprox(P.y, -1.8996, atol=1e-3)

    P, intervals = shubert_piyavskii(x->-sum(k*sin((k+1)*x + k) for k in 1 : 5), -10.0, 10.0, 70.0, 0.01, 0.05)
    @test isapprox(P.y, -12.031, atol=1e-3)

    a, b = bisection(x->x-1, 0.0, 100.0, 51.0)
    @test a == 0
    @test b == 50

    a, b = bisection(x->x-1, 0.0, 100.0, 0.1)
    @test abs(a-1) < 0.1
    @test abs(b-1) < 0.1
    @test b > a
    @test abs(b - a) < 0.1

    f = x->x-1
    a, b = bracket_sign_change(f, 50.0, 49.0)
    @test a < b
    @test f(a)*f(b) < 0
end