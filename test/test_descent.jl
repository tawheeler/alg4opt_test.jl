let
    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]
    ∇f = x -> (A + A')'*x
    H = x -> (A + A')'

    @test norm(line_search(f, Float64[ 0, 0], Float64[1,1]) - [0, 0]) ≤ 1e-10
    @test norm(line_search(f, Float64[-1,-1], Float64[1,1]) - [0, 0]) ≤ 1e-10
    @test norm(line_search(f, Float64[-1, 0], Float64[1,0]) - [0, 0]) ≤ 1e-10
    @test norm(line_search(f, Float64[-1,-1], Float64[1,0]) -
        (Float64[-1,-1] +  Optim.optimize(α->f(Float64[-1,-1] + Float64[α,0]), -2.0, 2.0).minimizer*Float64[1,0])) ≤ 1e-9

    x = [-2, -1.5]
    d = [1.0,1.0]
    β = 1e-4
    α = backtracking_line_search(f, ∇f, x, d, 100.0, 0.5, β)
    @test f(x + α*d) ≤ f(x) + β*α*dot(∇f(x), d)

    @test norm(trust_region_descent(f, ∇f, H, x, 4) - [0,0]) ≤ 1e-4
end