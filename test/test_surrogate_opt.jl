let
    X = [[1.0,1.0],
         [0.0,1.0],
         [1.0,0.0],
         [0.0,0.0]]

    @test μ(X, x->norm(x)) == [√2, 1, 1, 0]
    @test Σ(X, (x,x′)->norm(x-x′)) == [0.0 1.0 1.0  √2;
                                       1.0 0.0  √2 1.0;
                                       1.0  √2 0.0 1.0;
                                        √2 1.0 1.0 0.0]

    X′ = [[2.0,1.0],
          [1.0,1.0]]
    @test K(X, X′, (x,x′)->norm(x-x′)) == [1.0 0.0;
                                           2.0 1.0;
                                            √2 1.0;
                                            √5  √2]

    μ₁ = [1.0,2.0]
    Σ₁ = [1.0 -0.5; -0.5 2.0]
    srand(0)
    X = hcat(collect(mvnrand(μ₁, Σ₁) for i in 1 : 1000)...)
    D = fit(MvNormal, X)
    @test norm(D.μ - μ₁) ≤ 0.1
    @test norm(vec(full(D.Σ) - Σ₁)) ≤ 0.1

    GP = GaussianProcess(x->1.0, (x,x′)->exp(-norm(x-x′)), [[1.0,1.0]], [3.0], 0.5)
    rand(GP, X′)
    μₚ, νₚ = predict(GP, X′)
    @test all(μₚ .≈ [1.0, 1.0] + 4/3*[exp(-1), 1])
    @test all(νₚ .≈ [1.0, 1.0] - 2/3*[exp(-2), 1])

    @test prob_of_improvement( 0.0, 0.0, 1.0) ≈ 0.5
    @test prob_of_improvement(-1.0, 0.0, 1.0) ≈ cdf(Normal(0.0,1.0), -1.0)
    @test prob_of_improvement(-0.6,-0.5, 1.3) ≈ cdf(Normal(-0.5,1.3),-0.6)

    @test expected_improvement(0.0, 0.0, 1.0) ≈ pdf(Normal(0.0,1.0), 0.0)
    @test expected_improvement(-0.6,-0.5, 1.3) ≈ (-0.6 - (-0.5))*prob_of_improvement(-0.6,-0.5, 1.3) + 1.3*pdf(Normal(-0.5,1.3), -0.6)

    srand(0)
    X = [[1.0,1.0],
         [0.0,1.0],
         [1.0,0.0],
         [0.0,0.0]]
    GP = GaussianProcess(x->1.0, (x,x′)->exp(-norm(x-x′)), Vector{Float64}[], Float64[], 0.5)
    u_best, i_best = safe_opt(GP, X, 1, x->norm(x),   5.0)
    @test i_best == 4

    srand(0)
    GP = GaussianProcess(x->1.0, (x,x′)->exp(-norm(x-x′)), Vector{Float64}[], Float64[], 0.5)
    u_best, i_best = safe_opt(GP, X, 1, x->norm(x), -15.0)
    @test isnan(u_best)
    @test i_best == 0
end

let
    function f_(x::Vector{Float64})
        return norm(x)
    end
    f_(x1::Float64, x2::Float64) = f_([x1,x2])

    xdomain = ( -3, 3)
    ydomain = ( -3, 3)
    x_target_a =[-2.0, 1.0]
    x_target_b =[-2.0, 1.25]
    y_max = max(f_(x_target_a) + 1.0, f_(x_target_b) + 1.0)

    GP = GaussianProcess(x->y_max + 0.5, (x,x′)->exp(-0.5norm(x-x′)), Vector{Float64}[], Float64[], 0.01)
    β = 3.0
    n = 51
    X = Array{Vector{Float64}}(n*n)
    i = 0
    for x1 in linspace(xdomain..., n)
        for x2 in linspace(ydomain..., n)
            X[i+=1] = [x1,x2]
        end
    end

    m = length(X)
    u, ℓ = fill(Inf, m), fill(-Inf, m)
    S, M, E = falses(m), falses(m), falses(m)

    srand(0)

    i = indmin(norm(x_target_a-x,2) for x in X)
    push!(GP, X[i], f_(X[i]))

    i = indmin(norm(x_target_b-x,2) for x in X)
    push!(GP, X[i], f_(X[i]))

    for k in 2 : 15
        update_confidence_intervals!(GP, X, u, ℓ, β)
        S[:] = u .≤ y_max
        compute_sets!(GP, S, M, E, X, u, ℓ, y_max, β)
        i = get_new_query_point(M, E, u, ℓ)
        push!(GP, X[i], f_(X[i]))
    end
end

