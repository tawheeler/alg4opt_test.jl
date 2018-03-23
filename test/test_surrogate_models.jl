let
    X = [[1.0], [2.0]]

    @test design_matrix(X) ≈ [1.0 1.0; 1.0 2.0]
    y = [3.0,3.5]
    f = linear_regression(X, y)
    for i in 1 : length(X)
        @test f(X[i]) ≈ y[i]
    end
    @test f([0.0]) ≈ 2.5

    bases = polynomial_bases_1d(1, 3)
    @test [b([2.0]) for b in bases] ≈ [1, 2, 2^2, 2^3]

    bases = polynomial_bases(2, 2)
    @test [b([2.0,3.0]) for b in bases] ≈ [1, 2, 2^2, 3, 3*2, 3^2]

    f = x -> 0.5 + 2x[1] - 1(x[1])^2
    X = [[-1.0], [0.0], [2.0]]
    y = f.(X)

    bases = polynomial_bases(1, 2)
    @test [b([2.0]) for b in bases] ≈ [1, 2, 2^2]

    fhat = regression(X, y, polynomial_bases(1, 2))
    @test fhat.(X) ≈ y
    @test fhat([1.0]) ≈ f([1.0])

    bases = sinusoidal_bases_1d(1, 2, [0.0], Float64[2π])
    @test [b(Float64[π]) for b in bases] ≈ [0.5, 0.0, -1.0, 0.0, 1.0]

    # ks      powers
    # (0,0)   [0,0]
    # (1,0)   [1,0]
    # (2,0)   [1,0]
    # (3,0)   [2,0]
    # (4,0)   [2,0]
    # (0,1)   [0,1]
    # (1,1)   [1,1]
    # (2,1)   [1,1]
    # (0,2)   [0,1]
    # (1,2)   [1,1]
    # (2,2)   [1,1]
    # (0,3)   [0,2]
    # (0,4)   [0,2]
    bases = sinusoidal_bases(2, [0.0,2π], Float64[2π,4π])
    @test [b(Float64[π,2π]) for b in bases] ≈ [1/4, 0, -1/2, 0, 1/2, 0, 0, 0, 1/2, 0, -1, 0, 1/2]

    f = x -> sin(x[1]) - 0.2cos(x[1]) + 2
    X = [[0.0], [1.0], [2.0], [3.1], [4.5]]
    y = f.(X)
    fhat = regression(X, y, sinusoidal_bases(1, [0.0], Float64[2π]))
    @test fhat.(X) ≈ y
    @test fhat([0.5]) ≈ f([0.5])

    X = [[1.0], [2.0], [3.0]]
    bases = radial_bases(r->r^2, X)
    @test [b([1.5]) for b in bases] ≈ [0.5^2, 0.5^2, 1.5^2]
    y = [-1.0, 5.0, 6.0]
    fhat = regression(X, y, bases)
    @test fhat.(X) ≈ y

    fhat2 = regression(X, y, bases, 0.5)
    @test norm(fhat.(X) - y) < 1e-6

    f = x->1+2x[1]
    X = [[1.0], [2.0], [3.0]]
    y = f.(X)
    fit = linear_regression
    metric = (f, X, y) -> begin
        m = length(X)
        return sqrt(sum((f(X[i]) - y[i])^2 for i in m)/m)
    end
    @test isapprox(train_and_validate(X, y, TrainTest([1,2], [3]), fit, metric), 0.0, atol=1e-14)

    H = holdout_partition(2)
    @test (H.train == [1] && H.test == [2]) || (H.train == [2] && H.test == [1])
    srand(0)
    @test isapprox(random_subsampling(X, y, fit, metric), 0.0, atol=1e-14)

    srand(0)
    sets = k_fold_cross_validation_sets(5, 5)
    @test length(sets) == 5
    @test any(Set{Int}(S.train) == Set{Int}([1,2,3,4]) && S.test[1] == 5 for S in sets)
    @test any(Set{Int}(S.train) == Set{Int}([2,3,4,5]) && S.test[1] == 1 for S in sets)
    @test any(Set{Int}(S.train) == Set{Int}([3,4,5,1]) && S.test[1] == 2 for S in sets)
    @test any(Set{Int}(S.train) == Set{Int}([4,5,1,2]) && S.test[1] == 3 for S in sets)
    @test any(Set{Int}(S.train) == Set{Int}([5,1,2,3]) && S.test[1] == 4 for S in sets)

    srand(0)
    sets = k_fold_cross_validation_sets(4, 2)
    @test length(sets) == 2

    srand(0)
    sets = k_fold_cross_validation_sets(3, 2)
    cv_err = cross_validation_estimate(X, y, sets, fit, metric)
    @test cv_err > 0.1

    sets = bootstrap_sets(1, 2)
    @test length(sets) == 2
    @test sets[1].train == [1]
    @test sets[1].test == 1:1
    sets = bootstrap_sets(100,1)
    @test length(sets[1].train) == 100
    @test sets[1].test == 1:100

    srand(0)
    sets = bootstrap_sets(length(X), 10)
    boot_err = bootstrap_estimate(X, y, sets, fit, metric)
    @test boot_err < cv_err

    srand(0)
    sets = bootstrap_sets(length(X), 10)
    loo_boot_err = leave_one_out_bootstrap_estimate(X, y, sets, fit, metric)

    srand(0)
    sets = bootstrap_sets(length(X), 10)
    loo_boot_err = bootstrap_632_estimate(X, y, sets, fit, metric)
end