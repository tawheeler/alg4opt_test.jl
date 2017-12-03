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
    @test prob_of_improvement(-0.6,-0.5, 1.3) ≈ cdf(Normal(-0.5,sqrt(1.3)),-0.6)

    @test expected_improvement(0.0, 0.0, 1.0) ≈ pdf(Normal(0.0,1.0), 0.0)
    @test expected_improvement(-0.6,-0.5, 1.3) ≈ (-0.6 - (-0.5))*prob_of_improvement(-0.6,-0.5, 1.3) + sqrt(1.3)*pdf(Normal(-0.5,sqrt(1.3)), -0.6)

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
    function _flower(x; a=1, b=1, c=4)::Float64
        if isapprox(norm(x), 0.0)
            return 0.0
        end
        return a*norm(x) + b*sin(c*atan2(x[2], x[1]))
    end
    f_(x::Vector{Float64}) = _flower(x)
    f_(x1::Float64, x2::Float64) = f_([x1,x2])

    xdomain = ( -3, 3)
    ydomain = ( -3, 3)
    y_max = 2.0

    GP = GaussianProcess(x->y_max + 0.5, (x,x′)->exp(-norm(x-x′)), Vector{Float64}[], Float64[], 0.01)
    β = 3.0
    n = 51
    X = Array{Vector{Float64}}(n*n)
    i = 0
    for x1 in linspace(xdomain..., n)
        for x2 in linspace(ydomain..., n)
            X[i+=1] = [x1,x2]
        end
    end
    i = indmin(norm([-2,1]-x,2) for x in X)

    m = length(X)
    u, ℓ = fill(Inf, m), fill(-Inf, m)
    S, M, E = falses(m), falses(m), falses(m)

    srand(0)
    push!(GP, X[i], f_(X[i]) + randn()*GP.ν)
    update_confidence_intervals!(GP, X, u, ℓ, β)
    S[:] = u .≤ y_max
    best_val, i_best = findmin(u[S])
    i_best = findfirst(cumsum(S), i_best)

    for k in 2 : prod((4,4))
        compute_sets!(S, M, E, X, u, ℓ, y_max)
        i = get_new_query_point(M, E, u, ℓ)

        push!(GP, X[i], f_(X[i]))
        update_confidence_intervals!(GP, X, u, ℓ, β);
        S[:] = u .≤ y_max
        best_val, i_best = findmin(u[S])
        i_best = findfirst(cumsum(S), i_best)
    end
end