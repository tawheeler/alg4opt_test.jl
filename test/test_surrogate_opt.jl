let
    μ₁ = [1.0,2.0]
    Σ₁ = [1.0 -0.5; -0.5 2.0]
    srand(0)
    X = hcat(collect(mvnrand(μ₁, Σ₁) for i in 1 : 1000)...)
    D = fit(MvNormal, X)
    @test norm(D.μ - μ₁) ≤ 0.1
    @test norm(vec(full(D.Σ) - Σ₁)) ≤ 0.1

    # X_train
    # predict()
end