let
    M = 40
    k_max = 50
    f = x -> norm(x)
    srand(0)
    @test sum(genetic_algorithm(f, rand_population_binary(M, 5), k_max, TruncationSelection(5), SinglePointCrossover(), BitwiseMutation(0.1))) ≤ 1
    srand(0)
    @test sum(genetic_algorithm(f, rand_population_binary(M, 5), k_max, TournamentSelection(5), TwoPointCrossover(),    BitwiseMutation(0.1))) ≤ 1
    srand(0)
    @test sum(genetic_algorithm(f, rand_population_binary(M, 5), k_max, RouletteWheelSelection(), UniformCrossover(),   BitwiseMutation(0.1))) ≤ 1
    srand(0)
    @test f(genetic_algorithm(f, rand_population_uniform(M, [-2.0, -2.0], [2.0,2.0]), k_max, TruncationSelection(5), SinglePointCrossover(), GaussianMutation(0.01))) ≤ 0.1
    srand(0)
    @test f(genetic_algorithm(f, rand_population_normal(M, [0.0, 0.0], [1.0,1.0]), k_max, TournamentSelection(5), TwoPointCrossover(),    GaussianMutation(0.1))) ≤ 0.1
    srand(0)
    @test f(genetic_algorithm(f, rand_population_cauchy(M, [0.0, 0.0], [1.0,1.0]), k_max, RouletteWheelSelection(), UniformCrossover(),   GaussianMutation(0.1))) ≤ 0.1
    srand(0)
    @test f(genetic_algorithm(f, rand_population_normal(M, [0.0, 0.0], [1.0,1.0]), k_max, RouletteWheelSelection(), InterpolationCrossover(0.5),   GaussianMutation(0.1))) ≤ 0.1

    srand(0)
    population = [rand(Uniform(-10.0, 10.0), 2) for i in 1 : 100]
    @test f(rosenbrock(firefly(rosenbrock, population, 50))) < 0.01
end