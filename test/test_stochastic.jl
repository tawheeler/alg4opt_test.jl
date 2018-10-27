struct EvoStratParams
    μ::Vector{Float64}
    A::Matrix{Float64} # Σ = A'A
end
struct EvoStratGradient
    ∇μ::Vector{Float64}
    ∇A::Matrix{Float64}
end

function ∇logp(x::Vector{Float64}, θ::EvoStratParams)
    μ, A = θ.μ, θ.A
    Σ = A'*A
    Q = inv(Σ)
    ∇logpμ = Σ\(x - μ)
    ∇logpΣ = 0.5*Q*(x-μ)*(x-μ)'*Q - Q.*0.5
    ∇logpA = A*(∇logpΣ + ∇logpΣ')
    return EvoStratGradient(∇logpμ, ∇logpA)
end

Base.rand(θ::EvoStratParams) = rand(MvNormal(θ.μ, θ.A'*θ.A))

Base.:+(a::EvoStratGradient, b::EvoStratGradient) = EvoStratGradient(a.∇μ + b.∇μ, a.∇A + b.∇A)
Base.:*(x::Real, ∇::EvoStratGradient) = EvoStratGradient(∇.∇μ.*x, ∇.∇A.*x)
Base.:/(∇::EvoStratGradient, x::Real) = EvoStratGradient(∇.∇μ./x, ∇.∇A./x)
Base.:-(a::EvoStratParams, b::EvoStratGradient) = EvoStratParams(a.μ - b.∇μ, a.A - b.∇A)

let
    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]

    _rosenbrock(x; a=1, b=5) = (a-x[1])^2 + b*(x[2] - x[1]^2)^2

    seed!(0)
    @test f(mesh_adaptive_direct_search(f, [-2, -1.5], 0.0001)) < 0.01

    seed!(0)
    @test f(simulated_annealing(f, [-2, -1.5], MvNormal(Matrix(1.0I, 2, 2)), k->50/k, 250)) < 0.01
    seed!(0)
    @test norm(simulated_annealing(_rosenbrock, [-2,-1.5], MvNormal(Matrix(1.0I, 2, 2)), k->50/k, 250) - [1,1]) < 0.25

    seed!(0)
    @test f(adaptive_simulated_annealing(f,[-2, -1.5], ones(2), 50.0, 0.01)) < 1e-6
    seed!(0)
    @test norm(adaptive_simulated_annealing(_rosenbrock, [-2,-1.5], ones(2), 50.0, 0.01, ns=50) - [1,1]) < 1e-3

    seed!(0)
    P = MvNormal([-0.5,-1.5],[1.0,1.0])
    @test f(mean(cross_entropy_method(f, P, 10))) < 1e-5
    seed!(0)
    P = MvNormal([-0.5,-1.5], [5.0,5.0])
    @test norm(mean(cross_entropy_method(_rosenbrock, P, 15)) - [1,1]) < 1e-1

    seed!(0)
    θ = EvoStratParams([-0.5,-0.5], diagm(0=>[1.0,1.0]))
    θ = natural_evolution_strategies(f, θ, 30, α=0.25)
    P = MvNormal(θ.μ, θ.A'*θ.A)
    @test norm(params(P)[1]) < 1e-1
    @test norm(Matrix(params(P)[2])) < 1e-2
    # NOTE: I am convinced the alg is correct, and tested it manually on rosenbrock.

    seed!(0)
    @test norm(covariance_matrix_adaptation(f, [-0.5,-0.5], 10)) < 1e-1

    seed!(0)
    # NOTE: rosenbrock by itself has very large outliers in terms of y, so we are squashing it.
    @test norm(covariance_matrix_adaptation(x->sqrt(_rosenbrock(x)+1), [-2,-1.5], 20, m=20) - [1,1]) < 1e-1
end
let
    function _minimize(M::DescentMethod, f, ∇f, x, n, ε=0.001)
        init!(M, f, ∇f, x)
        for i in 1 : n
            x′ = step!(M, f, ∇f, x)
            if norm(x - x′) < ε
                break
            end
            x = x′
        end
        return x
    end

    A = Float64[1 -0.9; -0.9 1]
    f = x -> (x'*A*x)[1]
    ∇f = x -> (A + A')'*x
    x = Float64[-2, -1.5]

    for M in [
            NoisyDescent(GradientDescent(0.5), k->1/k, 0),
            NoisyDescent(ConjugateGradientDescent(NaN,NaN), k->0.001, 0),
        ]

        seed!(0)
        @test f(_minimize(M, f, ∇f, x, 50)) < 0.1
    end
end
# let
#     P = x -> [0.3,0.6,0.1][x]
#     T = x -> mod1(x+(rand() < 0.5 ? 1 : 2),3)
#     probT = (x,x′) -> 0.5

#     seed!(0)
#     counts = zeros(3)
#     for i in 1 : 1000
#         counts[metropolis_hastings(1, P, T, probT, burnin)] += 1
#     end
#     @test isapprox(counts[1] / sum(counts), 0.3, atol=0.1)
#     @test isapprox(counts[2] / sum(counts), 0.6, atol=0.1)

#     f = x -> -x
#     seed!(0)
#     x = markov_chain_monte_carlo(f, 1, T, probT, burnin=burnin)
#     @test x == 3
# end