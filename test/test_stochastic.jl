let
    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]

    rosenbrock(x; a=1, b=5) = (a-x[1])^2 + b*(x[2] - x[1]^2)^2

    srand(0)
    @test f(simulated_annealing(f, [-2, -1.5], MvNormal(eye(2)), k->50/k, 250)) < 0.01
    srand(0)
    @test norm(simulated_annealing(rosenbrock, [-2,-1.5], MvNormal(eye(2)), k->50/k, 250) - [1,1]) < 0.25

    srand(0)
    @test f(adaptive_simulated_annealing(f,[-2, -1.5], ones(2), 50.0, 0.01)) < 1e-6
    srand(0)
    @test norm(adaptive_simulated_annealing(rosenbrock, [-2,-1.5], ones(2), 50.0, 0.01, ns=50) - [1,1]) < 1e-3
end
# let
#     P = x -> [0.3,0.6,0.1][x]
#     T = x -> mod1(x+(rand() < 0.5 ? 1 : 2),3)
#     probT = (x,x′) -> 0.5

#     srand(0)
#     counts = zeros(3)
#     for i in 1 : 1000
#         counts[metropolis_hastings(1, P, T, probT, burnin)] += 1
#     end
#     @test isapprox(counts[1] / sum(counts), 0.3, atol=0.1)
#     @test isapprox(counts[2] / sum(counts), 0.6, atol=0.1)

#     f = x -> -x
#     srand(0)
#     x = markov_chain_monte_carlo(f, 1, T, probT, burnin=burnin)
#     @test x == 3
# end