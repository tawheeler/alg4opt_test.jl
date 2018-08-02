using Test
using Libdl
# using Distributions
# using ExprRules
using LightGraphs
# using SCS

import Base: rand
import Base.MathConstants: φ
import Statistics: var
import StatsBase: sample
import LinearAlgebra: ⋅, dot, norm, normalize, eye, I, diag, diagm, Diagonal, normalize!, triu, pinv
import Random: srand, randperm, bitrand
import Iterators: product
import Optim

# the -f option will cause fatal errors to error out runtests
fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"

# the -q option will quiet out error printing
quiet = length(ARGS) > 0 && ARGS[1] == "-q"

using DataFrames

function minimize(f::Function, a::Real, b::Real)
    a = convert(Float64, a)
    b = convert(Float64, b)
    return Optim.optimize(f, a, b).minimizer
end
function minimize(f::Function, x::Real)
    return Optim.optimize(f, x-100, x+100).minimizer
end
function minimize(f::Function, x::Vector{Float64})
    return Optim.optimize(f, x).minimizer
end

include(joinpath(@__DIR__, "..", "src", "all_julia_code.jl"))

function Base.push!(GP::GaussianProcess, x::Vector{Float64}, y::Float64)
    push!(GP.X, x)
    push!(GP.y, y)
    return GP
end
function Base.pop!(GP::GaussianProcess)
    pop!(GP.X)
    pop!(GP.y)
    return GP
end

# function minimize_lp_cp(LP)
#     A, b, c = LP.A, LP.b, LP.c
#     m, n = size(A)
#     z = ones(m)
#     Z = diagm([j ≥ 0? 1 : -1 for j in b])

#     A′ = hcat(A, Z)
#     b′ = b
#     c′ = vcat(zeros(n), z)
#     LP_init = LinearProgram(A′, b′, c′)
#     B = collect(1:m).+n
#     minimize_lp!(B, LP_init)

#     if any(i-> i > n, B)
#         error("infeasible")
#     end

#     A′′ = [A eye(m); zeros(m,n) eye(m)]
#     b′′ = vcat(b, zeros(m))
#     c′′ = c′
#     LP_opt = LinearProgram(A′′, b′′, c′′)
#     minimize_lp!(B, LP_opt)

#     x = get_vertex(B, LP_opt)[1:n]
#     b_inds = sort!(B)
#     n_inds = sort!(setdiff(1:n, B))

#     filter!(x->1 ≤ x ≤ n, b_inds)
#     filter!(x->1 ≤ x ≤ n, n_inds)

#     return x, b_inds, n_inds
# end

my_tests = [
    "test_derivatives.jl",
    "test_bracketing.jl",
    "test_descent.jl",
    "test_first_order.jl",
    "test_second_order.jl",
    "test_direct.jl",
    "test_stochastic.jl", #
    "test_population.jl",
    "test_penalty.jl",
    "test_linear.jl",
    "test_multiobjective.jl",
    "test_sampling_plans.jl",
    "test_surrogate_models.jl",
    "test_surrogate_opt.jl",
    # "test_design_under_uncertainty.jl",
    # "test_uncertaintyprop.jl",
    # "test_discrete.jl",
    # "test_expr.jl",
    # "test_mdo.jl",
    # "test_math.jl",
    # "test_test_functions.jl",
    ]

println("Running tests:")

anyerrors = false
for my_test in my_tests
    try
        include(my_test)
        println("\t\033[1m\033[32mPASSED\033[0m: $(my_test)")
    catch e
        global anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(stdout, e, backtrace())
            println()
        end
    end
end

if anyerrors
    throw("Tests failed")
end

println("DONE")