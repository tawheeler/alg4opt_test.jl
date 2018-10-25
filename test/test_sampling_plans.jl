let
    X = samples_full_factorial([-1.0, 1.0], [1.0, 2.0], [3,3])
    @test [-1.0,1.0] ∈ X
    @test [ 0.0,1.0] ∈ X
    @test [ 1.0,1.0] ∈ X
    @test [-1.0,1.5] ∈ X
    @test [ 0.0,1.5] ∈ X
    @test [ 1.0,1.5] ∈ X
    @test [-1.0,2.0] ∈ X
    @test [ 0.0,2.0] ∈ X
    @test [ 1.0,2.0] ∈ X
    @test length(X) == 9

    X = uniform_projection_plan(5, 3)

    X = Vector{Float64}[]
    push!(X, [1.0, 1.0])
    push!(X, [0.0, 1.0])
    push!(X, [1.0, 0.0])
    @test min_dist([2.0,2.0], X, 1) ≈ 2
    @test min_dist([2.0,2.0], X, 2) ≈ √(2)
    @test min_dist([2.0,2.0], X, Inf) ≈ 1
    @test min_dist([2.0,0.0], X, 2) ≈ 1

    X2 = Vector{Float64}[]
    push!(X2, [2.0,2.0])
    push!(X2, [3.0,3.0])
    @test d_max(X2, X, 2) ≈ hypot(2,2)


    X = Vector{Float64}[]
    for i in 1 : 100
        push!(X, randn(5))
    end
    @test length(greedy_local_search(X, 10)) == 10
    @test length(exchange_algorithm(X, 10)) == 10
    @test length(multistart_local_search(X, 10, greedy_local_search, 5)) == 10

    X = Vector{Float64}[]
    push!(X, [0,0])
    push!(X, [1,1])
    push!(X, [2,2])
    push!(X, [3,3])
    push!(X, [0,3])
    seed!(0)
    @test d_max(X, greedy_local_search(X, 2)) ≈ hypot(2,1)
    seed!(0)
    @test d_max(X, exchange_algorithm(X, 2)) ≈ hypot(2,1)

    X = Vector{Float64}[]
    push!(X, [0.0])
    push!(X, [1.0])
    push!(X, [5.0])
    push!(X, [6.0])
    dists = pairwise_distances(X)
    @test sort(dists) == sort(Float64[1,5,6,4,5,1])

    X = Vector{Float64}[]
    push!(X, [1.0,1.0])
    push!(X, [2.0,2.0])
    push!(X, [3.0,3.0])
    push!(X, [4.0,4.0])

    X2 = Vector{Float64}[]
    push!(X2, [1.0,3.0])
    push!(X2, [2.0,1.0])
    push!(X2, [3.0,2.0])
    push!(X2, [4.0,4.0])
    @test compare_sampling_plans(X, X2) ==  1
    @test compare_sampling_plans(X2, X) == -1
    @test compare_sampling_plans(X, X) == 0

    seed!(0)
    X3 = mutate!(deepcopy(X2))
    @test X2 != X3
    @test length(X2) == length(X3)

    X = get_filling_set_additive_recurrence(10)
    X = get_filling_set_additive_recurrence(10, 3)
    X = get_filling_set_halton(10)
    X = get_filling_set_halton(10,2)

    @test halton(1, 2) ≈ 1/2
    @test halton(2, 2) ≈ 1/4
    @test halton(3, 2) ≈ 3/4
    @test halton(4, 2) ≈ 1/8
    @test halton(5, 2) ≈ 5/8
    @test halton(6, 2) ≈ 3/8
    @test halton(7, 2) ≈ 7/8
    @test halton(8, 2) ≈ 1/16
    @test halton(9, 2) ≈ 9/16
    @test halton(1, 3) ≈ 1/3
    @test halton(2, 3) ≈ 2/3
    @test halton(3, 3) ≈ 1/9
    @test halton(4, 3) ≈ 4/9
    @test halton(5, 3) ≈ 7/9
    @test halton(6, 3) ≈ 2/9
    @test halton(7, 3) ≈ 5/9
    @test halton(8, 3) ≈ 8/9
    @test halton(9, 3) ≈ 1/27

    X = Vector{Float64}[]
    for i in 1 : 5
        for j in 1 : 5
            push!(X, [i,j])
        end
    end
    @test phiq(X, 3, 2) ≈ 3.9814903699475557
    @test phiq(X, 4, 3) ≈ 2.7670117482325867
end