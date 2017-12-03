using Base.Test
using Distributions
using ExprRules
using LightGraphs
using SCS

import Base: rand
import StatsBase: sample
import Optim

# the -f option will cause fatal errors to error out runtests
fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"

# the -q option will quiet out error printing
quiet = length(ARGS) > 0 && ARGS[1] == "-q"

using Base.Test
using DataFrames
using Compat

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

include(Pkg.dir("alg4opt_test", "src", "all_julia_code.jl"))

function Base.push!(GP::GaussianProcess, x::Vector{Float64}, y::Float64)
    push!(GP.X, x)
    push!(GP.y, y)
    return GP
end
function Base.pop!(GP::GaussianProcess, x::Vector{Float64}, y::Float64)
    pop!(GP.X)
    pop!(GP.y)
    return GP
end

my_tests = [
    "test_derivatives.jl",
    "test_bracketing.jl",
    "test_descent.jl",
    "test_first_order.jl",
    "test_second_order.jl",
    "test_direct.jl",
    "test_stochastic.jl",
    "test_population.jl",
    "test_penalty.jl",
    "test_linear.jl",
    "test_multiobjective.jl",
    "test_sampling_plans.jl",
    "test_surrogate_models.jl",
    "test_surrogate_opt.jl",
    "test_design_under_uncertainty.jl",
    "test_uncertaintyprop.jl",
    "test_discrete.jl",
    "test_expr.jl",
    "test_mdo.jl",
    "test_math.jl",
    ]

println("Running tests:")

anyerrors = false
for my_test in my_tests
    try
        include(my_test)
        println("\t\033[1m\033[32mPASSED\033[0m: $(my_test)")
    catch e
        anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(STDOUT, e, backtrace())
            println()
        end
    end
end

if anyerrors
    throw("Tests failed")
end

println("DONE")
