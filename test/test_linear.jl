let
    c = Float64[-5,-4,-3,0,0,0]
    A = Float64[2 3 1 1 0 0; 4 1 2 0 1 0; 3 4 2 0 0 1]
    b = Float64[5, 11, 8]
    LP = LinearProgram(A, b, c)
    β = [4,5,6]
    minimize_lp!(β, LP)
    x = get_vertex(β, LP)
    @test norm(x - [2,0,1,0,1,0]) < 0.01

    A = Float64[1 1 1 0; -4 2 0 1]
    b = Float64[9, 2]
    c = Float64[3, -1, 0, 0]
    LP = LinearProgram(A, b, c)
    β = [3,4]
    minimize_lp!(β, LP)
    x = get_vertex(β, LP)
    @test norm(x - [0,1,8,0]) < 0.0001

    A = Float64[1 1 -1; -1 2 0; 1 2 3]
    b = Float64[1, -2, 5]
    c = Float64[1, 1, -1]
    LP = LinearProgram(A, b, c)
    x = Float64[2, 0, 1]
    @test dual_certificate(LP, x, Float64[1, 0, 0], 1e-1)

    A = Float64[1   1 1 0;
                2 0.5 0 1]
    b = Float64[5, 8]
    c = Float64[-3, -2, 0, 0]
    LP = LinearProgram(A, b, c)
    x = minimize_lp(LP)
    @test x[3:4] == [0,0]
    @test isapprox(x[1], 11/3, atol=1e-8)
    @test isapprox(x[2],  4/3, atol=1e-8)
end