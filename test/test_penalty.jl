let
    @test norm(penalty_method(x->norm(x), x->(x[1]-2)^2, [0.0,10.0], 10) - [2,0]) < 0.01
    @test norm(augmented_lagrange_method(x->norm(x), x->[x[1] - 2, x[2] - 2], [0.0,0.0], 10) - [2,2]) < 0.01
end