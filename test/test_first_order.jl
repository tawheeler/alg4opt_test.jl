let
    function _minimize(M::DescentMethod, f, ∇f, x, n, ε=0.001)
        init!(M, f, ∇f, x)
        for i in 1 : n
            x′ = step(M, f, ∇f, x)
            if norm(x - x′) < ε
                break
            end
            x = x′
        end
        return x
    end

    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]
    ∇f = x -> (A + A')'*x
    x = Float64[-2, -1.5]

    for M in [
            GradientDescent(0.5),
            ConjugateGradientDescent(NaN,NaN),
            Momentum(0.1, 0.7, NaN),
            NesterovMomentum(0.1, 0.7, NaN),
            Adagrad(0.2, 1e-8, NaN),
            RMSprop(0.2, 0.45, 1e-8, NaN),
            Adadelta(0.2, 0.45, 0.45, 1e-8, NaN),
            Adam(0.2, 0.9, 0.9, 1e-8, 0, NaN, NaN),
            HyperGradientDescent(0.2, 1e-6, NaN, NaN),
            HyperNesterovMomentum(0.2, 1e-6, 0.93, NaN, NaN, NaN),
        ]
        
        @test f(_minimize(M, f, ∇f, x, 50)) < 0.1
    end
end