let
    f′ = directional_derivative(x->[2x[1], 2x[2]], [1.0,1.0], [1.0, 0.0])
    @test f′(0.0) ≈ 2.0
end