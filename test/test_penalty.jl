let
    @test norm(penalty_method(x->norm(x), x->max( 2-x[1],0)^2, [0.0,10.0], 10) - [2,0]) < 0.01
    @test norm(penalty_method(x->norm(x), x->max(-1-x[1],0)^2, [0.0,10.0], 10) - [0,0]) < 0.01
    @test norm(augmented_lagrange_method(x->norm(x), x->[x[1] - 2, x[2] - 2], [0.0,0.0], 10) - [2,2]) < 0.01
    @test norm(interior_point_method(x->norm(x), x->1/max(x[1] - 2,0), [3.0,1.0]) - [2,0]) < 0.01
    @test norm(interior_point_method(x->norm(x), x->1/max(x[1] + 1,0), [3.0,1.0]) - [0,0]) < 0.01
end