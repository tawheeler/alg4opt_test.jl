let
    function _minimize(M::DescentMethod, f, ∇f, x, n, ε=sqrt(eps()))
        x = deepcopy(x)
        init!(M, f, ∇f, x)
        for i in 1 : n
            @assert !any(x->isnan(x) || isinf(x), x)
            x′ = step!(M, f, ∇f, x)
            if norm(x - x′) < ε
                break
            end
            x = x′
        end
        return x
    end

    A = Float64[1 0; 0 1]
    f = x -> (x'*A*x)[1]
    ∇f = x -> (A + A')'*x
    H = x -> (A + A')'
    x = Float64[-2, -1.5]

    # Newton's should converge in one step
    @test norm(f(newtons_method(∇f, H, x, 0.5, 1))) ≈ 0
    @test f(_minimize(DFP(NaN), f, ∇f, x, 5)) < 0.001
    @test f(_minimize(BFGS(NaN), f, ∇f, x, 5)) < 0.001
    @test f(_minimize(LimitedMemoryBFGS(10,NaN,NaN,NaN), f, ∇f, x, 1)) < 0.001

    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]
    ∇f = x -> (A + A')'*x
    H = x -> (A + A')'
    x = Float64[-2, -1.5]

    # Newton's should converge in one step
    @test norm(f(newtons_method(∇f, H, x, 0.5, 1))) ≈ 0
    @test f(_minimize(DFP(NaN), f, ∇f, x, 5)) < 0.001
    @test f(_minimize(BFGS(NaN), f, ∇f, x, 5)) < 0.001
    @test f(_minimize(LimitedMemoryBFGS(10,NaN,NaN,NaN), f, ∇f, x, 5)) < 0.001
    @test f(_minimize(LimitedMemoryBFGS(2,NaN,NaN,NaN), f, ∇f, x, 5)) < 0.001

    f = x -> (1-x[1])^2 + 5*(4x[2] - x[1]^2)^2
    ∇f = x -> [2*(10x[1]^3 - 40x[1]*x[2] + x[1] - 1), -40*(x[1]^2 - 4x[2])]
    @test f(_minimize(DFP(NaN), f, ∇f, x, 15)) < 0.001
    @test f(_minimize(BFGS(NaN), f, ∇f, x, 15)) < 0.001
    @test f(_minimize(LimitedMemoryBFGS(5,NaN,NaN,NaN), f, ∇f, x, 10)) < 0.001


    f = x -> x^2
    f′ = x -> 2x
    @test f(secant_method(f′, 4, 3, 0.01)) < f(0.01)
end