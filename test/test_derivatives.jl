let

    for (f,x,∂) in [(x->x,    0.0, 1.0),
                    (x->x,    1.0, 1.0),
                    (x->x,    1.0, 1.0),
                    (x->x^2,  0.0, 0.0),
                    (x->x^2,  1.0, 2.0),
                    (x->x^2, -1.0,-2.0)]
        @test isapprox(derivative_forward(f, x), ∂, atol=1e-6)
        @test isapprox(derivative_central(f, x), ∂, atol=1e-6)
        @test isapprox(derivative_backward(f, x), ∂, atol=1e-6)
        @test isapprox(derivative_complex(f, x), ∂, atol=1e-6)
    end

    f′ = directional_derivative(x->[2x[1], 2x[2]], [1.0,1.0], [1.0, 0.0])
    @test f′(0.0) ≈ 2.0
end