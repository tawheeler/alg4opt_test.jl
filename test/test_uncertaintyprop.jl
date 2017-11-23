let
    x = [0.0]
    z = [1.0, 2.0]
    f = z -> sin(x[1]+z[1])*cos(x[1] + z[2])
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5])
    @test isapprox(ret[1], 0.0, atol=1e-10)
    @test isapprox(ret[2], 1.0, atol=1e-10)

    x = [-1.6]
    f = z -> sin(x[1]+z[1])*cos(x[1] + z[2])
    ret = taylor_approx(f, [0.0,0.0], [1.0,0.5])
    @test isapprox(ret[1], 0.0292, atol=1e-4)
    @test isapprox(ret[2], 0.4991, atol=1e-4)
end