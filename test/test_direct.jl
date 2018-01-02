let
    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]

    for x in [[-2,-1.5], [1.0,0.0], [0.0,0.0]]
        @test f(cyclic_coordinate_descent(f, x, 0.001)) < 0.01
        @test f(cyclic_coordinate_descent_with_acceleration_step(f, x, 0.001)) < 0.01
        @test f(powell(f, x, 0.001)) < 0.01
        @test f(generalized_pattern_search(f, x, 1.0, [[1.0,0.0],[-1.0,0.0],[0.0,1.0],[0.0,-1.0]], 0.001, 0.5)) < 0.01
        @test f(generalized_pattern_search(f, x, 1.0, [[1.0,0.0],[0.0,1.0],[-1.0,-1.0]], 0.001, 0.5)) < 0.01
        @test f(hooke_jeeves(f, x, 1.0, 0.001)) < 0.01
    end

    S = [
         [-2.0, -2.0],
         [-4.0, -4.0],
         [-4.0, -2.0],
        ]
    @test f(nelder_mead(f, S, 0.0001, α=1.0, β=2.0, γ=0.5)) < 0.01
    S = [
         [2.0, 2.0],
         [4.0, 4.0],
         [4.0, 2.0],
        ]
    @test f(nelder_mead(f, S, 0.001)) < 0.01

    @test f(direct(f, [-3.0, -3.0], [6.0, 7.0], 0.01, 10)) < 0.01
    @test abs(direct(x->sin(x[1]) + sin(2x[1]) + sin(4x[1]) + sin(8x[1]), [-2.0], [2.0], 0.001, 10)[1] - (-0.272)) < 0.001

    f2 = reparameterize_to_unit_hypercube(x->x[1], [2.0], [3.0])
    @test f2([0.0]) == 2.0
    @test f2([0.5]) == 2.5
    @test f2([1.0]) == 3.0
end

let
    # minimum close to [1,1]
    f = x -> exp(-dot(x,x)) - exp(-dot(x-[1,1],x-[1,1]))
    S = [
         [ 1.0, 1.0],
         [-1.0, 1.0],
         [ 0.0,-2.0],
        ]
    @test f(nelder_mead(f, S, 0.001)) < -0.85
end