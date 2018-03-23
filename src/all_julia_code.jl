#################### derivatives 1
function directional_derivative(∇f, x, d)
	α -> ∇f(x + α*d)⋅d
end
####################

#################### derivatives 2
diff_forward(f, x; h=sqrt(eps(Float64))) = (f(x+h) - f(x))/h
diff_central(f, x; h=cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/h
diff_backward(f, x; h=sqrt(eps(Float64))) = (f(x) - f(x-h))/h

####################

#################### derivatives 3
diff_complex(f, x; h=1e-20) = imag(f(x + h*im)) / h
####################

#################### bracketing 1
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end
####################

#################### bracketing 2
function fibonacci_search(f, a, b, n; ϵ=0.01)
    s = (1-√5)/(1+√5)
    ρ = 1 / (φ*(1-s^(n+1))/(1-s^n))
    d = ρ*b + (1 - ρ)*a
    yd = f(d)
    for i in 1 : n - 1
        if i == n - 1
            c = ϵ*a + (1-ϵ)*d
        else
            c = ρ*a + (1 - ρ)*b
        end
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
        ρ = 1 / (φ*(1-s^(n-i+1))/(1-s^(n-i)))
    end
    return a < b ? (a, b) : (b, a)
end
####################

#################### bracketing 3
function golden_section_search(f, a, b, n)
    ρ = φ-1
    d = ρ * b + (1 - ρ)*a
    yd = f(d)
    for i = 1 : n-1
        c = ρ*a + (1 - ρ)*b
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
    end
    return a < b ? (a, b) : (b, a)
end
####################

#################### bracketing 4
function quadratic_fit_search(f, a, b, c, n)
	ya, yb, yc = f(a), f(b), f(c)
	for i in 1:n-3
		x = 0.5*(ya*(b^2-c^2)+yb*(c^2-a^2)+yc*(a^2-b^2)) /
		        (ya*(b-c)    +yb*(c-a)    +yc*(a-b))
		yx = f(x)
		if x > b
			if yx > yb
				c, yc = x, yx
			else
				a, ya, b, yb = b, yb, x, yx
			end
		elseif x < b
			if yx > yb
				a, ya = x, yx
			else
				c, yc, b, yb = b, yb, x, yx
	        end
	    end
	end
	return (a, b, c)
end
####################

#################### bracketing 5
struct Pt
	x
	y
end
function _get_sp_intersection(A, B, l)
    t = ((A.y - B.y) - l*(A.x - B.x)) / 2l
    return Pt(A.x + t, A.y - t*l)
end
function shubert_piyavskii(f, a, b, l, ϵ, δ=0.01)

    m = (a+b)/2
    A, M, B = Pt(a, f(a)), Pt(m, f(m)), Pt(b, f(b))
    pts = [A, _get_sp_intersection(A, M, l),
    	   M, _get_sp_intersection(M, B, l), B]
    Δ = Inf
    while Δ > ϵ
		i = indmin(P.y for P in pts)
		P = Pt(pts[i].x, f(pts[i].x))
		Δ = P.y - pts[i].y

		P_prev = _get_sp_intersection(pts[i-1], P, l)
		P_next = _get_sp_intersection(P, pts[i+1], l)

		deleteat!(pts, i)
		insert!(pts, i, P_next)
		insert!(pts, i, P)
		insert!(pts, i, P_prev)
    end

    intervals = []
    i = 2*(indmin(P.y for P in pts[1:2:end])) - 1
    for j in 2:2:length(pts)
        if pts[j].y < pts[i].y
            dy = pts[i].y - pts[j].y
            x_lo = max(a, pts[j].x - dy/l)
            x_hi = min(b, pts[j].x + dy/l)
            if !isempty(intervals) && intervals[end][2] + δ ≥ x_lo
            	intervals[end] = (intervals[end][1], x_hi)
            else
            	push!(intervals, (x_lo, x_hi))
            end
        end
    end
    return (pts[i], intervals)
end
####################

#################### bracketing 6
function bisection(f′, a, b, ϵ)
    if a > b; a,b = b,a; end # ensure a < b

    ya, yb = f′(a), f′(b)
    if ya == 0; b = a; end
    if yb == 0; a = b; end

    while b - a > ϵ
        x = (a+b)/2
        y = f′(x)
        if y == 0
            a, b = x, x
        elseif sign(y) == sign(ya)
            a = x
        else
            b = x
        end
    end

    return (a,b)
end
####################

#################### bracketing 7
function bracket_sign_change(f′, a, b; k=2)
    if a > b; a,b = b,a; end # ensure a < b

    center, half_width = (b+a)/2, (b-a)/2
    while f′(a)*f′(b) > 0
        half_width *= k
        a = center - half_width
        b = center + half_width
    end

    return (a,b)
end
####################

#################### descent 1
function line_search(f, x, d)
    if norm(d) ≈ 0; return x; end; objective = α -> f(x + α*d)
    a, b = bracket_minimum(objective)
    α = minimize(objective, a, b)
    return x + α*d
end
####################

#################### descent 2
function backtracking_line_search(f, ∇f, x, d, α, p=0.5, β=1e-4)
	y, g = f(x), ∇f(x)
	while f(x + α*d) > y + β*α*(g⋅d)
		α *= p
	end
	α
end
####################

#################### descent 3
function strong_backtracking(f, ∇, x, d; α=1, β=1e-4, σ=0.1)
    y0, g0, y_prev, α_prev = f(x), ∇(x)⋅d, NaN, 0
    αlo, αhi = NaN, NaN

    # bracket phase
    while true
        y = f(x + α*d)
        if y > y0 + β*α*g0 || (!isnan(y_prev) && y ≥ y_prev)
            αlo, αhi = α_prev, α
            break
        end
        g = ∇(x + α*d)⋅d
        if abs(g) ≤ -σ*g0
            return α
        elseif g ≥ 0
            αlo, αhi = α, α_prev
            break
        end
        y_prev, α_prev, α = y, α, 2α
    end

    # zoom phase
    ylo = f(x + αlo*d)
    while true
        α = (αlo + αhi)/2
        y = f(x + α*d)
        if y > y0 + β*α*g0 || y ≥ ylo
            αhi = α
        else
            g = ∇(x + α*d)⋅d
            if abs(g) ≤ -σ*g0
                return α
            elseif g*(αhi - αlo) ≥ 0
                αhi = αlo
            end
            αlo = α
        end
    end
end
####################

#################### descent 4
function trust_region_descent(f, ∇f, H, x, k_max;
	η1=0.25, η2=0.5, γ1=0.5, γ2=2.0, δ=1.0)
	y = f(x)
	for k in 1 : k_max
		r = 0.0
		while r < η1
			x′, y′ = solve_trust_region_subproblem(∇f, H, x, δ)
			r = (y - f(x′)) / (y - y′)
			if r < η1
				δ *= γ1
			else
				x, y = x′, y′
				if r > η2
					δ *= γ2
				end
			end
		end
	end
	return x
end

using Convex
function solve_trust_region_subproblem(∇f, H, x0, δ)
	x = Variable(length(x0))
	p = Convex.minimize(∇f(x0)⋅(x-x0) + quadform(x-x0, H(x0))/2)
	p.constraints += norm(x-x0) <= δ
	solve!(p, SCSSolver(verbose=false), verbose=false)
	return (x.value, p.optval)
end
####################

#################### first-order 1
abstract type DescentMethod end
struct GradientDescent <: DescentMethod
	α
end
init!(M::GradientDescent, f, ∇f, x) = M
function step(M::GradientDescent, f, ∇f, x)
	α, g = M.α, ∇f(x)
	return x - α*g
end
####################

#################### first-order 2
mutable struct ConjugateGradientDescent <: DescentMethod
	d
	r
end
function init!(M::ConjugateGradientDescent, f, ∇f, x)
	M.d = M.r = -∇f(x)
	return M
end
function step(M::ConjugateGradientDescent, f, ∇f, x)
	d, r = M.d, M.r
    r′ = -∇f(x)
    β = max(0, dot(r′, r′-r)/(r⋅r))
    d′ = r′ + β*d
    x′ = line_search(f, x, d′)
    M.d, M.r = d′, r′
    return x′
end
####################

#################### first-order 3
mutable struct Momentum <: DescentMethod
	α # learning rate
	β # momentum decay
	v # momentum
end
function init!(M::Momentum, f, ∇f, x)
	M.v = zeros(length(x))
	return M
end
function step(M::Momentum, f, ∇f, x)
	α, β, v, g = M.α, M.β, M.v, ∇f(x)
	v[:] = β*v - α*g
	return x + v
end
####################

#################### first-order 4
mutable struct NesterovMomentum <: DescentMethod
	α # learning rate
	β # momentum decay
	v # momentum
end
function init!(M::NesterovMomentum, f, ∇f, x)
	M.v = zeros(length(x))
	return M
end
function step(M::NesterovMomentum, f, ∇f, x)
	α, β, v = M.α, M.β, M.v
	v[:] = β*v - α*∇f(x + β*v)
	return x + v
end
####################

#################### first-order 5
mutable struct Adagrad <: DescentMethod
	α # learning rate
	ϵ # small value
	s # sum of square gradient
end
function init!(M::Adagrad, f, ∇f, x)
	M.s = zeros(length(x))
	return M
end
function step(M::Adagrad, f, ∇f, x)
	α, ϵ, s, g = M.α, M.ϵ, M.s, ∇f(x)
	s[:] += g.*g
	return x - α*g ./ (sqrt.(s) + ϵ)
end
####################

#################### first-order 6
mutable struct RMSprop <: DescentMethod
	α # learning rate
	γ # decay
	ϵ # small value
	s # sum of square gradient
end
function init!(M::RMSprop, f, ∇f, x)
	M.s = zeros(length(x))
	return M
end
function step(M::RMSprop, f, ∇f, x)
	α, γ, ϵ, s, g = M.α, M.γ, M.ϵ, M.s, ∇f(x)
	s[:] = γ*s + (1-γ)*(g.*g)
	return x - α*g ./ (ϵ + sqrt.(s))
end
####################

#################### first-order 7
mutable struct Adadelta <: DescentMethod
	γs # gradient decay
	γx # update decay
	ϵ # small value
	s # sum of square gradients
	u # sum of square updates
end
function init!(M::Adadelta, f, ∇f, x)
	M.s = zeros(length(x))
	M.u = zeros(length(x))
	return M
end
function step(M::Adadelta, f, ∇f, x)
	γs, γx, ϵ, s, u, g = M.γs, M.γx, M.ϵ, M.s, M.u, ∇f(x)
	s[:] = γs*s + (1-γs)*g.*g
	Δx = - (ϵ + sqrt.(u)) ./ (ϵ + sqrt.(s)) .* g
	u[:] = γx*u + (1-γx)*Δx.*Δx
	return x + Δx
end
####################

#################### first-order 8
mutable struct Adam <: DescentMethod
	α # learning rate
	γv # decay
	γs # decay
	ϵ # small value
	k # step counter
	v # 1st moment estimate
	s # 2nd moment estimate
end
function init!(M::Adam, f, ∇f, x)
	M.k = 0
	M.v = zeros(length(x))
	M.s = zeros(length(x))
	return M
end
function step(M::Adam, f, ∇f, x)
	α, γv, γs, ϵ, k = M.α, M.γv, M.γs, M.ϵ, M.k
	s, v, g = M.s, M.v, ∇f(x)
	v[:] = γv*v + (1-γv)*g
	s[:] = γs*s + (1-γs)*g.*g
	M.k = k += 1
	v_hat = v ./ (1 - γv^k)
	s_hat = s ./ (1 - γs^k)
	return x - α*v_hat ./ (ϵ + sqrt.(s_hat))
end
####################

#################### first-order 9
mutable struct HyperGradientDescent <: DescentMethod
	α0 # initial learning rate
	μ # learning rate of the learning rate
	α # current learning rate
	g_prev # previous gradient
end
function init!(M::HyperGradientDescent, f, ∇f, x)
	M.α = M.α0
	M.g_prev = zeros(length(x))
	return M
end
function step(M::HyperGradientDescent, f, ∇f, x)
	α, μ, g, g_prev = M.α, M.μ, ∇f(x), M.g_prev
	α = α + μ*(g⋅g_prev)
	M.g_prev, M.α = g, α
	return x - α*g
end
####################

#################### first-order 10
mutable struct HyperNesterovMomentum <: DescentMethod
	α0 # initial learning rate
	μ # learning rate of the learning rate
	β # momentum decay
	v # momentum
	α # current learning rate
	g_prev # previous gradient
end
function init!(M::HyperNesterovMomentum, f, ∇f, x)
	M.α = M.α0
	M.v = zeros(length(x))
	M.g_prev = zeros(length(x))
	return M
end
function step(M::HyperNesterovMomentum, f, ∇f, x)
	α, β, μ = M.α, M.β, M.μ
	v, g, g_prev = M.v, ∇f(x), M.g_prev
	α = α - μ*(g⋅(-g_prev - β*v))
	v[:] = β*v + g
	M.g_prev, M.α = g, α
	return x - α*(g + β*v)
end
####################

#################### second-order 1
function newtons_method(∇f, H, x, ϵ, k_max)
	k, Δ = 1, Inf
	while norm(Δ) > ϵ && k ≤ k_max
		Δ = H(x) \ ∇f(x)
		x -= Δ
		k += 1
	end
	return x
end
####################

#################### second-order 2
function secant_method(f′, x0, x1, ϵ)
    g0 = f′(x0)
    Δ = Inf
    while abs(Δ) > ϵ
        g1 = f′(x1)
        Δ = (x1 - x0)/(g1 - g0)*g1
        x0, x1, g0 = x1, x1 - Δ, g1
    end
    return x1
end
####################

#################### second-order 3
mutable struct DFP <: DescentMethod
	Q
end
function init!(M::DFP, f, ∇f, x)
	M.Q = eye(length(x))
	return M
end
function step(M::DFP, f, ∇f, x)
	Q, g = M.Q, ∇f(x)
	x′ = line_search(f, x, -Q*g)
	g′ = ∇f(x′)
	δ = x′ - x
    γ = g′ - g
    Q[:] = Q - Q*γ*γ'*Q/(γ'*Q*γ) + δ*δ'/(δ'*γ)
    return x′
end
####################

#################### second-order 4
mutable struct BFGS <: DescentMethod
	Q
end
function init!(M::BFGS, f, ∇f, x)
	M.Q = eye(length(x))
	return M
end
function step(M::BFGS, f, ∇f, x)
	Q, g = M.Q, ∇f(x)
	x′ = line_search(f, x, -Q*g)
	g′ = ∇f(x′)
	δ = x′ - x
    γ = g′ - g
    Q[:] = Q - (δ*γ'*Q + Q*γ*δ')/(δ'*γ) +
               (1 + (γ'*Q*γ)/(δ'*γ))[1]*(δ*δ')/(δ'*γ)
    return x′
end
####################

#################### second-order 5
mutable struct LimitedMemoryBFGS <: DescentMethod
	m
	δs
	γs
end
function init!(M::LimitedMemoryBFGS, f, ∇f, x)
	M.δs = []
	M.γs = []
	return M
end
function step(M::LimitedMemoryBFGS, f, ∇f, x)
    δs, γs, g = M.δs, M.γs, ∇f(x)
    m = length(δs)
    q = g
    if m > 0
        for i in m : -1 : 1
            q -= (δs[i]⋅q)/(γs[i]⋅δs[i])*γs[i]
        end
        z = (γs[1] .* δs[1] .* q) / (γs[1]⋅γs[1])
        for i in 1 : m
            z += δs[i]*(δs[i]⋅q - γs[i]⋅z)/(γs[i]⋅δs[i])
        end
        x′ = line_search(f, x, -z)
    else
        x′ = line_search(f, x, -g)
    end
    g′ = ∇f(x′)
    push!(δs, x′ - x); push!(γs, g′ - g)
    while length(δs) > M.m
        shift!(δs); shift!(γs)
    end
    return x′
end
####################

#################### direct 1
basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]
####################

#################### direct 2
function cyclic_coordinate_descent(f, x, ϵ)
    Δ, n = Inf, length(x)
    while abs(Δ) > ϵ
    	x′ = copy(x)
    	for i in 1 : n
        	d = basis(i, n)
        	x = line_search(f, x, d)
        end
        Δ = norm(x - x′)
    end
    return x
end
####################

#################### direct 3
function cyclic_coordinate_descent_with_acceleration_step(f, x, ϵ)
    Δ, n = Inf, length(x)
    while abs(Δ) > ϵ
    	x′ = copy(x)
    	for i in 1 : n
        	d = basis(i, n)
        	x = line_search(f, x, d)
        end
        x = line_search(f, x, x - x′) # acceleration step
        Δ = norm(x - x′)
    end
    return x
end
####################

#################### direct 4
function powell(f, x, ϵ)
    n = length(x)
    U = [basis(i,n) for i in 1 : n]
    Δ = Inf
    while Δ > ϵ
        x′ = x
        for i in 1 : n
            d = U[i]
            x′ = line_search(f, x′, d)
        end
        for i in 1 : n-1
            U[i] = U[i+1]
        end
        U[n] = d = x′ - x
        x′ = line_search(f, x, d)
        Δ = norm(x′ - x)
        x = x′
    end
    return x
end
####################

#################### direct 5
function hooke_jeeves(f, x, α, ϵ, γ=0.5)
    y, n = f(x), length(x)
    while α > ϵ
        improved = false
        for i in 1 : n
            x′ = x + α*basis(i, n)
            y′ = f(x′)
            if y′ < y
                x, y, improved = x′, y′, true
            else
                x′ = x - α*basis(i, n)
                y′ = f(x′)
                if y′ < y
                    x, y, improved = x′, y′, true
                end
            end
        end

        if !improved
            α *= γ
        end
    end
    return x
end
####################

#################### direct 6
function generalized_pattern_search(f, x, α, D, ϵ, γ=0.5)
    y, n = f(x), length(x)
    while α > ϵ
    	improved = false
        for (i,d) in enumerate(D)
            x′ = x + α*d
            y′ = f(x′)
            if y′ < y
                x, y, improved = x′, y′, true
                D = unshift!(deleteat!(D, i), d)
                break
            end
        end
        if !improved
            α *= γ
        end
    end
    return x
end
####################

#################### direct 7
function nelder_mead(f, S, ϵ; α=1.0, β=2.0, γ=0.5)
    Δ, y_arr = Inf, f.(S)
    while Δ > ϵ
        p = sortperm(y_arr) # sort lowest to highest
        S, y_arr = S[p], y_arr[p]
        xl, yl = S[1], y_arr[1] # lowest
        xh, yh = S[end], y_arr[end] # highest
        xs, ys = S[end-1], y_arr[end-1] # second-highest
        xm = mean(S[1:end-1]) # centroid
        xr = xm + α*(xm - xh) # reflection point
        yr = f(xr)

        if yr < yl
            xe = xm + β*(xr-xm) # expansion point
            ye = f(xe)
            S[end],y_arr[end] = ye < yr ? (xe, ye) : (xr, yr)
        elseif yr > ys
            if yr ≤ yh
                xh, yh, S[end], y_arr[end] = xr, yr, xr, yr
            end
            xc = xm + γ*(xh - xm) # contraction point
            yc = f(xc)
            if yc > yh
                for i in 2 : length(y_arr)
                    S[i] = (S[i] + xl)/2
                    y_arr[i] = f(S[i])
                end
            else
                S[end], y_arr[end] = xc, yc
            end
        else
            S[end], y_arr[end] = xr, yr
        end

        Δ = std(y_arr, corrected=false)
    end
    return S[indmin(y_arr)]
end
####################

#################### direct 8
function direct(f, a, b, ϵ, k_max)
    g = reparameterize_to_unit_hypercube(f, a, b)
    intervals = Intervals()
    n = length(a)
    c = fill(0.5, n)
    interval = Interval(c, g(c), fill(0, n))
    add_interval!(intervals, interval)
    c_best, y_best = copy(interval.c), interval.y

    for k in 1 : k_max
        S = get_opt_intervals(intervals, ϵ, y_best)
        to_add = Interval[]
        for interval in S
            append!(to_add, divide(g, interval))
            dequeue!(intervals[min_depth(interval)])
        end
        for interval in to_add
            add_interval!(intervals, interval)
            if interval.y < y_best
                c_best, y_best = copy(interval.c), interval.y
            end
        end
    end

    return rev_unit_hypercube_parameterization(c_best, a, b)
end
####################

#################### direct 9
rev_unit_hypercube_parameterization(x, a, b) = x.*(b-a) + a
function reparameterize_to_unit_hypercube(f, a, b)
    Δ = b-a
    return x->f(x.*Δ + a)
end
####################

#################### direct 10
using DataStructures
struct Interval
    c
    y
    depths
end
min_depth(interval) = minimum(interval.depths)
const Intervals = Dict{Int,PriorityQueue{Interval, Float64}}
function add_interval!(intervals, interval)
	d = min_depth(interval)
    if !haskey(intervals, d)
        intervals[d] = PriorityQueue{Interval, Float64}()
    end
    return enqueue!(intervals[d], interval, interval.y)
end
####################

#################### direct 11
function get_opt_intervals(intervals, ϵ, y_best)
    max_depth = maximum(keys(intervals))
    stack = [DataStructures.peek(intervals[max_depth])[1]]
    d = max_depth-1
    while d ≥ 0
        if haskey(intervals, d) && !isempty(intervals[d])
            interval = DataStructures.peek(intervals[d])[1]
            x, y = 0.5*3.0^(-min_depth(interval)), interval.y

            while !isempty(stack)
            	interval1 = stack[end]
            	x1 = 0.5*3.0^(-min_depth(interval1))
            	y1 = interval1.y
            	l1 = (y - y1)/(x - x1)
            	if y1 - l1*x1 > y_best - ϵ || y < y1
            		pop!(stack)
            	elseif length(stack) > 1
            		interval2 = stack[end-1]
            		x2 = 0.5*3.0^(-min_depth(interval2))
            		y2 = interval2.y
            		l2 = (y1 - y2)/(x1 - x2)
            		if l2 > l1
            			pop!(stack)
                    else
                        break
            		end
                else
                    break
            	end
            end

            push!(stack, interval) # add new point
        end
        d -= 1
    end
    return stack
end
####################

#################### direct 12
function divide(f, interval)
    c, d, n = interval.c, min_depth(interval), length(interval.c)
    dirs = find(interval.depths .== d)
    cs = [(c + 3.0^(-d-1)*basis(i,n),
           c - 3.0^(-d-1)*basis(i,n)) for i in dirs]
    vs = [(f(C[1]), f(C[2])) for C in cs]
    minvals = [min(V[1], V[2]) for V in vs]

    intervals = Interval[]
    depths = copy(interval.depths)
    for j in sortperm(minvals)
        depths[dirs[j]] += 1
        C, V = cs[j], vs[j]
        push!(intervals, Interval(C[1], V[1], copy(depths)))
        push!(intervals, Interval(C[2], V[2], copy(depths)))
    end
    push!(intervals, Interval(c, interval.y, copy(depths)))
    return intervals
end
####################

#################### stochastic 1
mutable struct NoisyDescent <: DescentMethod
	submethod
	σ
	k
end
function init!(M::NoisyDescent, f, ∇f, x)
	init!(M.submethod, f, ∇f, x)
	M.k = 1
	return M
end
function step(M::NoisyDescent, f, ∇f, x)
	x = step(M.submethod, f, ∇f, x)
	σ = M.σ(M.k)
	x += σ.*randn(length(x))
	M.k += 1
	return x
end
####################

#################### stochastic 2
function rand_positive_spanning_set(α, n)
    δ = round(Int, 1/sqrt(α))
    L = diagm(δ*rand([1,-1], n))
    for i in 1 : n-1
    	for j in i+1:n
    		L[i,j] = rand(-δ+1:δ-1)
    	end
    end
    D = L[randperm(n),:]
    D = L[:,randperm(n)]
    D = hcat(D, -sum(D,2))
    return [D[:,i] for i in 1 : n+1]
end
####################

#################### stochastic 3
function mesh_adaptive_direct_search(f, x, ϵ)
    α, y, n = 1, f(x), length(x)
    while α > ϵ
    	improved = false
        for (i,d) in enumerate(rand_positive_spanning_set(α, n))
            x′ = x + α*d
            y′ = f(x′)
            if y′ < y
                x, y, improved = x′, y′, true
				x′ = x + 3α*d
				y′ = f(x′)
				if y′ < y
					x, y = x′, y′
				end
                break
            end
        end
        α = improved ? min(4α, 1) : α/4
    end
    return x
end
####################

#################### stochastic 4
function simulated_annealing(f, x, T, t, k_max)
    y = f(x)
    x_best, y_best = x, y
    for k in 1 : k_max
        x′ = x + rand(T)
        y′ = f(x′)
        Δy = y′ - y
        if Δy ≤ 0 || rand() < exp(-Δy/t(k))
            x, y = x′, y′
        end
        if y′ < y_best
            x_best, y_best = x′, y′
        end
    end
    return x_best
end
####################

#################### stochastic 5
function corana_update!(v, a, c, ns)
    for i in 1 : length(v)
        ai, ci = a[i], c[i]
        if ai > 0.6ns
        	v[i] *= (1 + ci*(ai/ns - 0.6)/0.4)
        elseif ai < 0.4ns
        	v[i] /= (1 + ci*(0.4-ai/ns)/0.4)
        end
    end
    return v
end
####################

#################### stochastic 6
function adaptive_simulated_annealing(f, x, v, t, ϵ;
    ns=20, nϵ=4, nt=max(100,5length(x)),
    γ=0.85, c=fill(2,length(x)) )

    y = f(x)
    x_best, y_best = x, y
    y_arr, n, U = [], length(x), Uniform(-1.0,1.0)
    a,counts_cycles,counts_resets = zeros(n), 0, 0

    while true
	    for i in 1:n
	        x′ = x + basis(i,n)*rand(U)*v[i]
	        y′ = f(x′)
	        Δy = y′ - y
	        if Δy < 0 || rand() < exp(-Δy/t)
	            x, y = x′, y′
	            a[i] += 1
	            if y′ < y_best; x_best, y_best = x′, y′; end
	        end
	    end

	    counts_cycles += 1
	    counts_cycles ≥ ns || continue

	    counts_cycles = 0
	    corana_update!(v, a, c, ns)
	    fill!(a, 0)
	    counts_resets += 1
	    counts_resets ≥ nt || continue

	    t *= γ
	    counts_resets = 0
	    push!(y_arr, y)

	    if !(length(y_arr) > nϵ && y_arr[end] - y_best ≤ ϵ &&
	         all(abs(y_arr[end]-y_arr[end-u]) ≤ ϵ for u in 1:nϵ))
	        x, y = x_best, y_best
	    else
	    	break
	    end
	end
    return x_best
end
####################

#################### stochastic 7
using Distributions
function cross_entropy_method(f, P, k_max, m=100, m_elite=10)
    for k in 1 : k_max
        samples = rand(P, m)
        order = sortperm([f(samples[:,i]) for i in 1:m])
        P = fit(typeof(P), samples[:,order[1:m_elite]])
    end
    return P
end
####################

#################### stochastic 8
using Distributions
function natural_evolution_strategies(f, θ, k_max; m=100, α=0.01)
    for k in 1 : k_max
        population = [rand(θ) for i in 1 : m]
        θ -= α*sum(f(x)*∇logp(x, θ) for x in population)/m
    end
    return θ
end
####################

#################### stochastic 9
function covariance_matrix_adaptation(f, x, k_max;
	σ = 1.0,
	m = 4 + floor(Int, 3*log(length(x))),
	m_elite = div(m,2))

	μ, n = copy(x), length(x)
	ws = normalize!(vcat(log((m+1)/2) - [log(i) for i in 1:m_elite],
	                zeros(m - m_elite)), 1)
	μ_eff = 1 / sum(ws.^2)
	cσ = (μ_eff + 2)/(n + μ_eff + 5)
	dσ = 1 + 2max(0, sqrt((μ_eff-1)/(n+1))-1) + cσ
	cΣ = (4 + μ_eff/n)/(n + 4 + 2μ_eff/n)
	c1 = 2/((n+1.3)^2 + μ_eff)
	cμ = min(1-c1, 2*(μ_eff-2+1/μ_eff)/((n+2)^2 + μ_eff))
	E = n^0.5*(1-1/(4n)+1/(21*n^2))
	pσ, pΣ, Σ = zeros(n), zeros(n), eye(n)
	for k in 1 : k_max
	    P = MvNormal(μ, σ^2*Σ)
	    xs = [rand(P) for i in 1 : m]
	    ys = [f(x) for x in xs]
	    is = sortperm(ys) # best to worst

	    # selection and mean update
	    δs = [(x - μ)/σ for x in xs]
	    δw = sum(ws[i]*δs[is[i]] for i in 1 : m_elite)
	    μ += σ*δw

	    # step-size control
	    C = Σ^-0.5
	    pσ = (1-cσ)*pσ + sqrt(cσ*(2-cσ)*μ_eff)*C*δw
	    σ *= exp(cσ/dσ * (norm(pσ)/E - 1))

	    # covariance adaptation
	    hσ = norm(pσ)/sqrt(1-(1-cσ)^(2k)) < (1.4+2/(n+1))*E ? 1 : 0
	    pΣ = (1-cΣ)*pΣ + hσ*sqrt(cΣ*(2-cΣ)*μ_eff)*δw
	    w0 = [ws[i]≥0 ? ws[i] : n*ws[i]/norm(C*δs[is[i]])^2
	    	  for i in 1:m]
	    Σ = (1-c1-cμ) * Σ +
	        c1*(pΣ*pΣ' + (1-hσ) * cΣ*(2-cΣ) * Σ) +
	        cμ*sum(w0[i]*δs[is[i]]*δs[is[i]]' for i in 1 : m)
	    Σ = triu(Σ)+triu(Σ,1)' # enforce symmetry
	end
	return μ
end
####################

#################### population 1
function rand_population_uniform(m, a, b)
    d = length(a)
    return [a+rand(d).*(b-a) for i in 1:m]
end
####################

#################### population 2
using Distributions
function rand_population_normal(m, μ, Σ)
    D = MvNormal(μ,Σ)
    return [rand(D) for i in 1:m]
end
####################

#################### population 3
using Distributions
function rand_population_cauchy(m, μ, σ)
    n = length(μ)
    return [[rand(Cauchy(μ[j],σ[j])) for j in 1:n] for i in 1:m]
end
####################

#################### population 4
function genetic_algorithm(f, population, k_max, S, C, M)
    for k in 1 : k_max
        parents = select(S, f.(population))
        children = [crossover(C,population[p[1]],population[p[2]])
                    for p in parents]
        population .= mutate.(M, children)
    end
    population[indmin(f.(population))]
end
####################

#################### population 5
rand_population_binary(m, n) = [bitrand(n) for i in 1:m]
####################

#################### population 6
abstract type SelectionMethod end

struct TruncationSelection <: SelectionMethod
    k::Int # top k to keep
end
function select(t::TruncationSelection, y)
    p = sortperm(y)
    return [p[rand(1:t.k, 2)] for i in y]
end

struct TournamentSelection <: SelectionMethod
    k::Int
end
function select(t::TournamentSelection, y)
    getparent() = begin
        p = randperm(length(y))
        p[indmin(y[p[1:t.k]])]
    end
    return [[getparent(), getparent()] for i in y]
end

struct RouletteWheelSelection <: SelectionMethod end
function select(::RouletteWheelSelection, y)
    y = maximum(y) - y
    cat = Categorical(normalize(y, 1))
    return [rand(cat, 2) for i in y]
end
####################

#################### population 7
abstract type CrossoverMethod end
struct SinglePointCrossover <: CrossoverMethod end
function crossover(::SinglePointCrossover, a, b)
    i = rand(1:length(a))
    return vcat(a[1:i], b[i+1:end])
end

struct TwoPointCrossover <: CrossoverMethod end
function crossover(::TwoPointCrossover, a, b)
    n = length(a)
    i, j = rand(1:n, 2)
    if i > j
        (i,j) = (j,i)
    end
    return vcat(a[1:i], b[i+1:j], a[j+1:n])
end

struct UniformCrossover <: CrossoverMethod end
function crossover(::UniformCrossover, a, b)
    child = copy(a)
    for i in 1 : length(a)
        if rand() < 0.5
            child[i] = b[i]
        end
    end
    return child
end
####################

#################### population 8
struct InterpolationCrossover <: CrossoverMethod
    λ
end
crossover(C::InterpolationCrossover, a, b) = (1-C.λ)*a + C.λ*b
####################

#################### population 9
abstract type MutationMethod end
struct BitwiseMutation <: MutationMethod
	λ
end
function mutate(M::BitwiseMutation, child)
    return [rand() < M.λ ? !v : v for v in child]
end

struct GaussianMutation <: MutationMethod
    σ
end
function mutate(M::GaussianMutation, child)
    return child + randn(length(child))*M.σ
end
####################

#################### population 10
using StatsBase
function differential_evolution(f, population, k_max; p=0.5, w=1)
    n, m = length(population[1]), length(population)
    for k in 1 : k_max
        for (k,x) in enumerate(population)
            weights = Weights([j!=k for j in 1 : m])
            a, b, c = sample(population, weights, 3, replace=false)
            z = a + w*(b-c)
            j = rand(1:n)
            x′ = [i == j || rand() < p ? z[i] : x[i] for i in 1:n]
            if f(x′) < f(x)
                x[:] = x′
            end
        end
    end
    return population[indmin(f.(population))]
end
####################

#################### population 11
mutable struct Particle
    x
    v
    x_best
end
####################

#################### population 12
function particle_swarm_optimization(f, population, k_max;
    w=1, c1=1, c2=1)
    n = length(population[1].x)
    x_best, y_best = copy(population[1].x_best), Inf
    for P in population
        y = f(P.x)
        if y < y_best; x_best[:], y_best = P.x, y; end
    end
    for k in 1 : k_max
        for P in population
            r1, r2 = rand(n), rand(n)
            P.x += P.v
            P.v = w*P.v + c1*r1.*(P.x_best - P.x) +
                          c2*r2.*(x_best - P.x)
            y = f(P.x)
            if y < y_best; x_best[:], y_best = P.x, y; end
            if y < f(P.x_best); P.x_best[:] = P.x; end
        end
    end
    return population
end
####################

#################### population 13
using Distributions
function firefly(f, population, k_max; β=1, α=0.1, I=r->exp(-r^2))
    N = MvNormal(eye(length(population[1])))
    for k in 1 : k_max
        for a in population, b in population
            if f(b) < f(a)
                r = norm(b-a)
                a[:] += β*I(r)*(b-a) + α*rand(N)
            end
        end
    end
    return population[indmin([f(x) for x in population])]
end
####################

#################### population 14
using Distributions
mutable struct Nest
	x # position
	y # value, f(x)
end
function cuckoo_search(f, population, k_max; p_a=0.1, C=Cauchy(0,1))
	m, n = length(population), length(population[1].x)
    a = round(Int, m*p_a)
	for k in 1 : k_max
		i, j = rand(1:m), rand(1:m)
        x = population[j].x + [rand(C) for k in 1 : n]
        y = f(x)
        if y < population[i].y
            population[i].x[:] = x
            population[i].y = y
        end

        p = sortperm(population, by=nest->nest.y, rev=true)
        for i in 1 : a
            j = rand(1:m-a)+a
            population[p[i]] = Nest(population[p[j]].x +
                                 [rand(C) for k in 1 : n],
                                 f(population[p[i]].x)
                                 )
        end
	end
	return population
end
####################

#################### penalty 1
function penalty_method(f, p, x, k_max; ρ=1, γ=2)
	for k in 1 : k_max
		x = minimize(x -> f(x) + ρ*p(x), x)
		ρ *= γ
		if p(x) == 0
			return x
		end
	end
	return x
end
####################

#################### penalty 2
function augmented_lagrange_method(f, h, x, k_max; ρ=1, γ=2)
	λ = zeros(length(h(x)))
	for k in 1 : k_max
		p = x -> f(x) + ρ/2*sum(h(x).^2) - λ⋅h(x)
		x = minimize(x -> f(x) + p(x), x)
		ρ *= γ
		λ -= ρ*h(x)
	end
	return x
end
####################

#################### penalty 3
function interior_point_method(f, p, x; ρ=1, γ=2, ϵ=0.001)
	delta = Inf
	while delta > ϵ
		x′ = minimize(x -> f(x) + p(x)/ρ, x)
		delta = norm(x′ - x)
		x = x′
		ρ *= γ
	end
	return x
end
####################

#################### linear 1
mutable struct LinearProgram
    A
    b
    c
end
function get_vertex(B, LP)
    A, b, c = LP.A, LP.b, LP.c
    b_inds = sort!(collect(B))
    AB = A[:,b_inds]
    xB = AB\b
    x = zeros(length(c))
    x[b_inds] = xB
	return x
end
####################

#################### linear 2
function edge_transition(LP, B, q)
	A, b, c = LP.A, LP.b, LP.c
    n = size(A, 2)
    b_inds = sort(B)
    n_inds = sort!(setdiff(1:n, B))
    AB = A[:,b_inds]
    d, xB = AB\A[:,n_inds[q]], AB\b

	p, xq′ = 0, Inf
	for i in 1 : length(d)
	    if d[i] > 0
	        v = xB[i] / d[i]
	        if v < xq′
	            p, xq′ = i, v
	        end
	    end
	end
	return (p, xq′)
end
####################

#################### linear 3
function step_lp!(B, LP)
    A, b, c = LP.A, LP.b, LP.c
    n = size(A, 2)
    b_inds = sort!(B)
    n_inds = sort!(setdiff(1:n, B))
    AB, AV = A[:,b_inds], A[:,n_inds]
    xB = AB\b
    cB = c[b_inds]
    λ = AB' \ cB
    cV = c[n_inds]
    μV = cV - AV'*λ

    q, p, xq′, Δ = 0, 0, Inf, Inf
    for i in 1 : length(μV)
    	if μV[i] < 0
    		pi, xi′ = edge_transition(LP, B, i)
    		if μV[i]*xi′ < Δ
    			q, p, xq′, Δ = i, pi, xi′, μV[i]*xi′
    		end
    	end
   	end
   	if q == 0
        return (B, true) # optimal point found
    end

    if isinf(xq′)
        error("unbounded")
    end
    B[findfirst(B, b_inds[p])] = n_inds[q] # swap indices
    return (B, false) # new vertex but not optimal
end
####################

#################### linear 4
function minimize_lp!(B, LP)
    done = false
    while !done
        B, done = step_lp!(B, LP)
    end
    return B
end
####################

#################### linear 5
function minimize_lp(LP)
    A, b, c = LP.A, LP.b, LP.c
    m, n = size(A)
    z = ones(m)
    Z = diagm([j ≥ 0? 1 : -1 for j in b])

    A′ = hcat(A, Z)
    b′ = b
    c′ = vcat(zeros(n), z)
    LP_init = LinearProgram(A′, b′, c′)
    B = collect(1:m).+n
    minimize_lp!(B, LP_init)

	if any(i-> i > n, B)
		error("infeasible")
	end

	A′′ = [A eye(m); zeros(m,n) eye(m)]
	b′′ = vcat(b, zeros(m))
	c′′ = c′
	LP_opt = LinearProgram(A′′, b′′, c′′)
	minimize_lp!(B, LP_opt)
	return get_vertex(B, LP_opt)[1:n]
end
####################

#################### linear 6
function dual_certificate(LP, x, μ, ϵ=1e-6)
	A, b, c = LP.A, LP.b, LP.c
	primal_feasible = all(x .≥ 0) && A*x ≈ b
	dual_feasible = all(A'*μ .≤ c)
	return primal_feasible && dual_feasible &&
	       isapprox(c⋅x, b⋅μ, atol=ϵ)
end
####################

#################### multiobjective 1
dominates(y, y′) = all(y′ - y .≥ 0) && any(y′ - y .> 0)
####################

#################### multiobjective 2
function naive_pareto(xs, ys)
    pareto_xs, pareto_ys = similar(xs, 0), similar(ys, 0)
    for (x,y) in zip(xs,ys)
        if !any(dominates(y′,y) for y′ in ys)
            push!(pareto_xs, x)
            push!(pareto_ys, y)
        end
    end
    return (pareto_xs, pareto_ys)
end
####################

#################### multiobjective 3
function weight_pareto(f1, f2, npts)
    return [
        Optim.optimize(x->w1*f1(x) + (1-w1)*f2(x), 0, π/2).minimizer
        for w1 in linspace(0,1,npts)
    ]
end
####################

#################### multiobjective 4
function vector_evaluated_genetic_algorithm(f, population,
	k_max, S, C, M)
	m = length(f(population[1]))
	m_pop = length(population)
	m_subpop = m_pop ÷ m
    for k in 1 : k_max
    	ys = f.(population)
    	parents = select(S, [y[1] for y in ys])[1:m_subpop]
    	for i in 2 : m
    		subpop=select(S,[y[i] for y in ys])[1:m_subpop]
    		append!(parents, subpop)
    	end

    	p = randperm(2m_pop)
    	p_ind=i->parents[mod(p[i]-1,m_pop)+1][(p[i]-1)÷m_pop + 1]
    	parents = [[p_ind(i), p_ind(i+1)] for i in 1 : 2 : 2m_pop]
        children = [crossover(C,population[p[1]],population[p[2]])
                    for p in parents]
        population = [mutate(M, c) for c in children]
    end
    return population
end
####################

#################### multiobjective 5
function get_non_domination_levels(ys)
    L, m = 0, length(ys)
    levels = zeros(Int, m)
    while minimum(levels) == 0
        L += 1
        for (i,y) in enumerate(ys)
            if levels[i] == 0 &&
               !any((levels[i] == 0 || levels[i] == L) &&
                    dominates(ys[i],y) for i in 1 : m)
                levels[i] = L
            end
        end
    end
    return levels
end
####################

#################### multiobjective 6
function discard_closest_pair!(xs, ys)
	index, min_dist = 0, Inf
	for (i,y) in enumerate(ys)
		for (j, y′) in enumerate(ys[i+1:end])
			dist = norm(y - y′)
			if dist < min_dist
				index, min_dist = rand([i,j]), dist
			end
		end
	end
	deleteat!(xs, index)
	deleteat!(ys, index)
	return (xs, ys)
end
####################

#################### multiobjective 7
function update_pareto_filter!(filter_xs, filter_ys, xs, ys;
	capacity=length(xs),
	)
    for (x,y) in zip(xs, ys)
    	if !any(dominates(y′,y) for y′ in filter_ys)
            push!(filter_xs, x)
            push!(filter_ys, y)
        end
    end
    filter_xs, filter_ys = naive_pareto(filter_xs, filter_ys)
    while length(filter_xs) > capacity
    	discard_closest_pair!(filter_xs, filter_ys)
    end
    return (filter_xs, filter_ys)
end
####################

#################### sampling-plans 1
import IterTools: product
function samples_full_factorial(a, b, m)
	ranges = [linspace(a[i], b[i], m[i]) for i in 1 : length(a)]
    collect.(collect(product(ranges...)))
end
####################

#################### sampling-plans 2
function uniform_projection_plan(m, n)
	perms = [randperm(m) for i in 1 : n]
	[[perms[i][j] for i in 1 : n] for j in 1 : m]
end
####################

#################### sampling-plans 3
function pairwise_distances(X, p=2)
    m = length(X)
    [norm(X[i]-X[j], p) for i in 1:(m-1) for j in (i+1):m]
end
####################

#################### sampling-plans 4
function compare_sampling_plans(A, B, p=2)
	pA = sort(pairwise_distances(A, p))
	pB = sort(pairwise_distances(B, p))
	for (dA, dB) in zip(pA, pB)
		if dA < dB
			return 1
		elseif dA > dB
			return -1
		end
	end
	return 0
end
####################

#################### sampling-plans 5
function mutate!(X)
    m, n = length(X), length(X[1])
    j = rand(1:n)
    i = randperm(m)[1:2]
    X[i[1]][j], X[i[2]][j] = X[i[2]][j], X[i[1]][j]
    return X
end
####################

#################### sampling-plans 6
function phiq(X, q=1, p=2)
	dists = pairwise_distances(X, p)
	return sum(dists.^(-q))^(1/q)
end
####################

#################### sampling-plans 7
min_dist(a, B, p) = minimum(norm(a-b, p) for b in B)
d_max(A, B, p=2) = maximum(min_dist(a, B, p) for a in A)
####################

#################### sampling-plans 8
function greedy_local_search(X, m, d=d_max)
	S = [X[rand(1:m)]]
	for i in 2 : m
		j = indmin(x ∈ S ? Inf : d(X, push!(copy(S), x))
		           for x in X)
		push!(S, X[j])
	end
	return S
end
####################

#################### sampling-plans 9
function exchange_algorithm(X, m, d=d_max)
	S = X[randperm(m)]
	δ, done = d(X, S), false
	while !done
		best_pair = (0,0)
		for i in 1 : m
			s = S[i]
			for (j,x) in enumerate(X)
				if !in(x, S)
					S[i] = x
					δ′ = d(X, S)
					if δ′ < δ
						δ = δ′
						best_pair = (i,j)
					end
				end
			end
			S[i] = s
		end
		done = best_pair == (0,0)
		if !done
			i,j = best_pair
			S[i] = X[j]
		end
	end
	return S
end
####################

#################### sampling-plans 10
function multistart_local_search(X, m, alg, k_max, d=d_max)
	sets = [alg(X, m, d) for i in 1 : k_max]
	return sets[indmin(d(X, S) for S in sets)]
end
####################

#################### sampling-plans 11
using Primes
function get_filling_set_additive_recurrence(m; c=φ-1)
    X = [rand()]
    for i in 2 : m
        push!(X, mod(X[end] + c, 1))
    end
    return X
end
function get_filling_set_additive_recurrence(m, n)
    ps = primes(max(ceil(Int, n*(log(n) + log(log(n)))), 6))
    seqs = [get_filling_set_additive_recurrence(m, c=sqrt(p))
            for p in ps[1:n]]
    return [collect(x) for x in zip(seqs...)]
end
####################

#################### sampling-plans 12
using Primes
function halton(i, b)
    result, f = 0.0, 1.0
    while i > 0
        f = f / b;
        result = result + f * mod(i, b)
        i = floor(Int, i / b)
    end
    return result
end
get_filling_set_halton(m; b=2) = [halton(i,b) for i in 1: m]
function get_filling_set_halton(m, n)
    bs = primes(max(ceil(Int, n*(log(n) + log(log(n)))), 6))
    seqs = [get_filling_set_halton(m, b=b) for b in bs[1:n]]
    return [collect(x) for x in zip(seqs...)]
end
####################

#################### surrogate-models 1
function design_matrix(X)
	n, m = length(X[1]), length(X)
	return [j==0 ? 1.0 : X[i][j] for i in 1:m, j in 0:n]
end
function linear_regression(X, y)
	θ = pinv(design_matrix(X))*y
	return x -> θ⋅[1; x]
end
####################

#################### surrogate-models 2
function regression(X, y, bases)
    B = [b(x) for x in X, b in bases]
    θ = pinv(B)*y
    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
end
####################

#################### surrogate-models 3
import IterTools: product
polynomial_bases_1d(i, k) = [x->x[i]^p for p in 0:k]
function polynomial_bases(n, k)
	bases = [polynomial_bases_1d(i, k) for i in 1 : n]
	terms = Function[]
	for ks in product([0:k for i in 1:n]...)
		if sum(ks) ≤ k
			push!(terms,
				x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
		end
	end
	return terms
end
####################

#################### surrogate-models 4
function sinusoidal_bases_1d(j, k, a, b)
	T = b[j] - a[j]
	bases = Function[x->1/2]
	for i in 1 : k
		push!(bases, x->sin(2π*i*x[j]/T))
		push!(bases, x->cos(2π*i*x[j]/T))
	end
	return bases
end
function sinusoidal_bases(k, a, b)
	n = length(a)
	bases = [sinusoidal_bases_1d(i, k, a, b) for i in 1 : n]
	terms = Function[]
	for ks in product([0:2k for i in 1:n]...)
		powers = [div(k+1,2) for k in ks]
		if sum(powers) ≤ k
			push!(terms,
				x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
		end
	end
	return terms
end
####################

#################### surrogate-models 5
radial_bases(ψ, C, p=2) = [x->ψ(norm(x - c, p)) for c in C]
####################

#################### surrogate-models 6
function regression(X, y, bases, λ)
    B = [b(x) for x in X, b in bases]
    θ = (B'B + λ*I)\B'y
    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
end
####################

#################### surrogate-models 7
struct TrainTest
    train
    test
end
function train_and_validate(X, y, tt, fit, metric)
    model = fit(X[tt.train], y[tt.train])
    return metric(model, X[tt.test], y[tt.test])
end
####################

#################### surrogate-models 8
function holdout_partition(m, h=div(m,2))
    p = randperm(m)
    train = p[(h+1):m]
    holdout = p[1:h]
    return TrainTest(train, holdout)
end
####################

#################### surrogate-models 9
function random_subsampling(X, y, fit, metric;
	h=div(length(X),2), k_max=10)
	m = length(X)
    mean(train_and_validate(X, y, holdout_partition(m, h),
          fit, metric) for k in 1 : k_max)
end
####################

#################### surrogate-models 10
function k_fold_cross_validation_sets(m, k)
    perm = randperm(m)
    sets = TrainTest[]
    for i = 1:k
        validate = perm[i:k:m];
        train = perm[setdiff(1:m, i:k:m)]
        push!(sets, TrainTest(train, validate))
    end
    return sets
end
function cross_validation_estimate(X, y, sets, fit, metric)
	mean(train_and_validate(X, y, tt, fit, metric)
		  for tt in sets)
end
####################

#################### surrogate-models 11
bootstrap_sets(m, b) = [TrainTest(rand(1:m, m), 1:m) for i in 1 : b]
####################

#################### surrogate-models 12
function bootstrap_estimate(X, y, sets, fit, metric)
	mean(train_and_validate(X, y, tt, fit, metric) for tt in sets)
end
####################

#################### surrogate-models 13
function leave_one_out_bootstrap_estimate(X, y, sets, fit, metric)
	m, b = length(X), length(sets)
	ε = 0.0
	models = [fit(X[tt.train], y[tt.train]) for tt in sets]
	for j in 1 : m
		c = 0
		δ = 0.0
		for i in 1 : b
			if j ∉ sets[i].train
				c += 1
				δ += metric(models[i], [X[j]], [y[j]])
			end
		end
		ε += δ/c
	end
	return ε/m
end
####################

#################### surrogate-models 14
function bootstrap_632_estimate(X, y, sets, fit, metric)
	models = [fit(X[tt.train], y[tt.train]) for tt in sets]
	ϵ_loob = leave_one_out_bootstrap_estimate(X,y,sets,fit,metric)
	ϵ_boot = bootstrap_estimate(X,y,sets,fit,metric)
    return 0.632ϵ_loob + 0.368ϵ_boot
end
####################

#################### prob-surrogate-models 1
μ(X, m) = [m(x) for x in X]
Σ(X, k) = [k(x,x′) for x in X, x′ in X]
K(X, X′, k) = [k(x,x′) for x in X, x′ in X′]
####################

#################### prob-surrogate-models 2
mutable struct GaussianProcess
	m # mean
	k # covariance function
	X # design points
	y # objective values
	ν # noise variance
end
####################

#################### prob-surrogate-models 3
function mvnrand(μ, Σ, inflation=1e-6)
	N = MvNormal(μ, Σ + inflation*I)
	return rand(N)
end
Base.rand(GP, X) = mvnrand(μ(X, GP.m), Σ(X, GP.k))
####################

#################### prob-surrogate-models 4
function predict(GP, X_pred)
    m, k, ν = GP.m, GP.k, GP.ν
    tmp = K(X_pred, GP.X, k) / (K(GP.X, GP.X, k) + ν*I)
    μₚ = μ(X_pred, m) + tmp*(GP.y - μ(GP.X, m))
    S = K(X_pred, X_pred, k) - tmp*K(GP.X, X_pred, k)
    νₚ = diag(S) .+ eps() # eps prevents numerical issues
    return (μₚ, νₚ)
end
####################

#################### surrogate-optimization 1
prob_of_improvement(y_min, μ, σ) = cdf(Normal(μ, σ), y_min)
####################

#################### surrogate-optimization 2
function expected_improvement(y_min, μ, σ)
    p_imp = prob_of_improvement(y_min, μ, σ)
    p_ymin = pdf(Normal(μ, σ), y_min)
    return (y_min - μ)*p_imp + σ*p_ymin
end
####################

#################### surrogate-optimization 3
function safe_opt(GP, X, i, f, y_max; β=3.0, k_max=10)
    push!(GP, X[i], f(X[i])) # make first observation

    m = length(X)
    u, l = fill(Inf, m), fill(-Inf, m)
    S, M, E = falses(m), falses(m), falses(m)

    for k in 1 : k_max
        update_confidence_intervals!(GP, X, u, l, β)
        compute_sets!(S, M, E, X, u, l, y_max)
        i = get_new_query_point(M, E, u, l)
        i != 0 || break
        push!(GP, X[i], f(X[i]))
    end

    # return the best point
    update_confidence_intervals!(GP, X, u, l, β)
    S[:] = u .≤ y_max
    if any(S)
        u_best, i_best = findmin(u[S])
        i_best = findfirst(cumsum(S), i_best)
        return (u_best, i_best)
    else
        return (NaN,0)
    end
end
####################

#################### surrogate-optimization 4
function update_confidence_intervals!(GP, X, u, l, β)
    μₚ, νₚ = predict(GP, X)
    u[:] = μₚ + sqrt.(β*νₚ)
    l[:] = μₚ - sqrt.(β*νₚ)
    return (u, l)
end
####################

#################### surrogate-optimization 5
function compute_sets!(S, M, E, X, u, l, y_max)
	fill!(M, false)
    fill!(E, false)

    # safe set
    S[:] = u .≤ y_max

    if any(S)

        # potential maximizers
        M[S] = u[S] .≥ maximum(l[S])

        # maximum width (in M)
        w_max = maximum(u[M] - l[M])

        # expanders - skip values in M or those with w ≤ w_max
        E[:] = S .& .~M # skip points in M
        if any(E)
            E[E] = maximum(u[E] - l[E]) .> w_max
            for (i,e) in enumerate(E)
                if e && u[i] - l[i] > w_max
                    push!(GP, X[i], l[i])
                    μₚ, νₚ = predict(GP, X[.~S])
                    pop!(GP)
                    E[i] = any(μₚ + sqrt.(β*νₚ) ≥ y_max)
                    if E[i]; w_max = u[i] - l[i]; end
                end
            end
        end
    end

    return (S,M,E)
end
####################

#################### surrogate-optimization 6
function get_new_query_point(M, E, u, l)
    ME = M .| E
    if any(ME)
    	return findfirst(cumsum(ME), indmax(u[ME] - l[ME]))
    else
    	return 0
    end
end
####################

#################### uncertaintyprop 1
function taylor_approx(f, μ, ν, secondorder=false)
    μhat = f(μ)
    ∇ = (z -> ForwardDiff.gradient(f, z))(μ)
    νhat = ∇.^2⋅ν
    if secondorder
        H = (z -> ForwardDiff.hessian(f, z))(μ)
        μhat += (diag(H)⋅ν)/2
        νhat += ν⋅(H.^2*ν)/2
    end
    return (μhat, νhat)
end
####################

#################### uncertaintyprop 2
using Polynomials
function legendre(i)
	n = i-1
    p = Poly([-1,0,1])^n
    for i in 1 : n
        p = polyder(p)
    end
    return p / (2^n * factorial(n))
end
function laguerre(i)
    p = Poly([1])
    for j in 2 : i
        p = polyint(polyder(p) - p) + 1
    end
    return p
end
function hermite(i)
    p = Poly([1])
    x = Poly([0,1])
    for j in 2 : i
        p = x*p - polyder(p)
    end
    return p
end
####################

#################### uncertaintyprop 3
using Polynomials
function orthogonal_recurrence(bs, p, dom, ϵ=1e-6)
    i = length(bs)
    c1 = quadgk(z->z*bs[i](z)^2*p(z), dom..., abstol=ϵ)[1]
    c2 = quadgk(z->  bs[i](z)^2*p(z), dom..., abstol=ϵ)[1]
    α = c1 / c2
    if i > 1
        c3 = quadgk(z->bs[i-1](z)^2*p(z), dom..., abstol=ϵ)[1]
        β = c2 / c3
        return Poly([-α, 1])*bs[i] - β*bs[i-1]
    else
        return Poly([-α, 1])*bs[i]
    end
end
####################

#################### uncertaintyprop 4
function polynomial_chaos_bases(bases1d)
	bases = []
	for a in product(bases1d...)
		push!(bases,
			z -> prod(b(z[i]) for (i,b) in enumerate(a))
		)
	end
	return bases
end
####################

#################### uncertaintyprop 5
function bayesian_monte_carlo(GP, w, μz, Σz)
	W = diagm(w.^2)
	invK = inv(K(GP.X, GP.X, GP.k))
	q = [exp(-0.5*((z-μz)'inv(W+Σz)*(z-μz))[1]) for z in GP.X]
	q .*= (det(W\Σz + I))^(-0.5)
	μ = q'*invK*GP.y
	ν = (det(2W\Σz + I))^(-0.5) - (q'*invK*q)[1]
	return (μ, ν)
end
####################

#################### discrete 1
mutable struct IntegerLinearProgram
    A
    b
    c
end
relax(IP) = LinearProgram(IP.A, IP.b, IP.c)
round_ip(IP) = round.(Int, minimize_lp(relax(IP)))
####################

#################### discrete 2
isint(x, ϵ=1e-10) = abs(round(x) - x) ≤ ϵ
function is_totally_unimodular(A::Matrix)
    # all entries must be in [0,1,-1]
    if any(a ∉ (0,-1,1) for a in A)
        return false
    end
    # brute force check every subdeterminant
    r,c = size(A)
    for i in 1 : min(r,c)
        for a in IterTools.subsets(1:r, i)
            for b in IterTools.subsets(1:c, i)
                B = A[a,b]
                if det(B) ∉ (0,-1,1)
                    return false
                end
            end
        end
    end
    return true
end
function is_totally_unimodular(IP)
    return is_totally_unimodular(IP.A) &&
           all(isint, IP.b) && all(isint, IP.c)
end
####################

#################### discrete 3
frac(x) = modf(x)[1]
function cutting_plane(IP)
    LP = relax(IP)
    x, b_inds, v_inds = minimize_lp_cp(LP)
    n_orig = length(x)
    while !all(isint.(x))

        AB, AV = LP.A[:,b_inds], LP.A[:,v_inds]
        Abar = AB\AV

        b, n = 0, length(x)
        for i in 1 : n
            if !isint(x[i])
                b += 1
                A2 = [LP.A zeros(size(LP.A,1));
                      zeros(1,size(LP.A,2)+1)]
                A2[end,end] = 1
                A2[end,v_inds] = (x->floor(x) - x).(Abar[b,:])
                b2 = vcat(LP.b, -frac(x[i]))
                c2 = vcat(LP.c, 0)
                LP = LinearProgram(A2,b2,c2)
            end
        end
        x, b_inds, v_inds = minimize_lp_cp(LP)
    end

    return round.(Int, x[1:n_orig])
end
####################

#################### discrete 4
import DataStructures: PriorityQueue
function minimize_lp_and_y(LP)
    try
        x = minimize_lp(LP)
        return (x, x⋅LP.c)
    catch
        return (fill(NaN, length(LP.c)), Inf)
    end
end
function branch_and_bound(IP)
    LP = relax(IP)
    x, y = minimize_lp_and_y(LP)
    n = length(x)
    x_best, y_best = deepcopy(x), Inf
    Q = PriorityQueue()
    enqueue!(Q, (LP,x,y), y)
    while !isempty(Q)
        LP, x, y = dequeue!(Q)
        if any(isnan.(x)) || all(isint.(x[1:n]))
            if y < y_best
                x_best, y_best = x[1:n], y
            end
        else
            i = indmax(abs(v - round(v)) for v in x)
            # x_i ≤ floor(x_i)
            A, b, c = LP.A, LP.b, LP.c
            A2=[A zeros(size(A,1)); [j==i for j in 1:size(A,2)]' 1]
            b2, c2 = vcat(b, floor(x[i])), vcat(c, 0)
            LP2 = LinearProgram(A2,b2,c2)
            x2, y2 = minimize_lp_and_y(LP2)
            if y2 ≤ y_best
                enqueue!(Q, (LP2,x2,y2), y2)
            end
            # x_i ≥ ceil(x_i)
            A2=[A zeros(size(A,1)); [j==i for j in 1:size(A,2)]' -1]
            b2, c2 = vcat(b, ceil(x[i])), vcat(c, 0)
            LP2 = LinearProgram(A2,b2,c2)
            x2, y2 = minimize_lp_and_y(LP2)
            if y2 ≤ y_best
                enqueue!(Q, (LP2,x2,y2), y2)
            end
        end
    end
    return round.(Int, x_best)
end
####################

#################### discrete 5
function padovan_topdown(n, P=Dict())
    if !haskey(P, n)
        P[n] = n < 3 ? 1 :
               padovan_topdown(n-2,P) + padovan_topdown(n-3,P)
    end
    return P[n]
end
function padovan_bottomup(n)
    P = Dict(0=>1,1=>1,2=>1)
    for i in 3 : n
        P[i] = P[i-2] + P[i-3]
    end
    return P[n]
end
####################

#################### discrete 6
function knapsack(v, w, w_max)
    n = length(v)
    y = Dict((0,j) => 0.0 for j in 0:w_max)
    for i in 1 : n
        for j in 0 : w_max
            y[i,j] = w[i] > j ? y[i-1,j] :
                     max(y[i-1,j],
                         y[i-1,j-w[i]] + v[i])
        end
    end

    # recover solution
    x, j = falses(n), w_max
    for i in n: -1 : 1
        if w[i] ≤ j && y[i,j] - y[i-1, j-w[i]] == v[i]
            # the ith element is in the knapsack
            x[i] = true
            j -= w[i]
        end
    end
    return x
end
####################

#################### discrete 7
function edge_attractiveness(graph, τ, η; α=1, β=5)
    A = Dict()
    for i in 1 : nv(graph)
        neighbors = out_neighbors(graph, i)
        for j in neighbors
            v = τ[(i,j)]^α * η[(i,j)]^β
            A[(i,j)] = v
        end
    end
    return A
end
####################

#################### discrete 8
import StatsBase: Weights, sample
function run_ant(graph, lengths, τ, A, x_best, y_best)
    x = [1]
    while length(x) < nv(graph)
        i = x[end]
        neighbors = setdiff(out_neighbors(graph, i), x)
        if isempty(neighbors) # ant got stuck
            return (x_best, y_best)
        end

        as = [A[(i,j)] for j in neighbors]
        push!(x, neighbors[sample(Weights(as))])
    end

    l = sum(lengths[(x[i-1],x[i])] for i in 2:length(x))
    for i in 2 : length(x)
        τ[(x[i-1],x[i])] += 1/l
    end
    if l < y_best
        return (x, l)
    else
        return (x_best, y_best)
    end
end
####################

#################### discrete 9
function ant_colony_optimization(graph, lengths;
    m=1000, k_max=100, α=1.0, β=5.0, ρ=0.5,
    η=Dict((e.src,e.dst)=>1/lengths[(e.src,e.dst)]
            for e in edges(graph)))

    τ = Dict((e.src,e.dst)=>1.0 for e in edges(graph))
    x_best, y_best = [], Inf
    for k in 1 : k_max
    	A = edge_attractiveness(graph, τ, η, α=α, β=β)
        for (e,v) in τ
            τ[e] = (1-ρ)*v
        end
        for ant in 1 : m
            x_best,y_best = run_ant(graph,lengths,τ,A,x_best,y_best)
        end
    end
    return x_best
end
####################

#################### expr 1
function _nested_monte_carlo_expression_discovery(
    f, grammar, root, stack, n, xb, yb)

    if n < 1 || isempty(stack)
        for l in stack
            typ = return_type(grammar, l.parent)
            insert!(root, l, rand(RuleNode, grammar, typ))
        end
        y = f(root)
        if y < yb
            xb, yb = deepcopy(root), y
        end
    else
        l = pop!(stack)
        for i in grammar[child_types(grammar, l.parent)[l.i]]
            ctypes = child_types(grammar, i)
            if !isempty(ctypes)
                node2 = RuleNode(i, [RuleNode(0) for c in ctypes])
                insert!(root, l, node2)
                for i in length(ctypes) : -1 : 1
                    push!(stack, NodeLoc(node2, i))
                end
                xb, yb = _nested_monte_carlo_expression_discovery(
                            f, grammar, root, stack, n-1, xb, yb)
                for i in 1 : length(node2.children)
                    pop!(stack)
                end
            else
                node2 = RuleNode(i)
                xb, yb = _nested_monte_carlo_expression_discovery(
                            f, grammar, root, stack, n-1, xb, yb)
            end
        end
        push!(stack, l)
    end
    return xb, yb
end
####################

#################### expr 2
function nested_monte_carlo_expression_discovery(f, grammar, typ, n)
    xb, yb = RuleNode(0), Inf
    for i in grammar[typ]
        ctypes = child_types(grammar, i)
        if !isempty(ctypes)
            root = RuleNode(i, [RuleNode(0) for c in ctypes])
            stack = [NodeLoc(root,i) for i in length(ctypes):-1:1]
            xb, yb = _nested_monte_carlo_expression_discovery(
                            f, grammar, root, stack, n-1, xb, yb)
        else
            y = f(RuleNode(i))
            if y < yb
                xb, yb = deepcopy(RuleNode(i)), y
            end
        end
    end
    return xb
end
####################

#################### expr 3
struct TreeCrossover <: CrossoverMethod
    grammar
    max_depth
end
function crossover(C::TreeCrossover, a, b)
    child = deepcopy(a)
    crosspoint = sample(b)
    typ = return_type(C.grammar, crosspoint.ind)
    d_subtree = depth(crosspoint)
    d_max = C.max_depth + 1 - d_subtree
    if d_max > 0 && contains_returntype(child, C.grammar, typ, d_max)
        loc = sample(NodeLoc, child, typ, C.grammar, d_max)
        insert!(child, loc, deepcopy(crosspoint))
    end
    child
end
####################

#################### expr 4
struct TreeMutation <: MutationMethod
    grammar
    p
end
function mutate(M::TreeMutation, a)
    child = deepcopy(a)
    if rand() < M.p
        loc = sample(NodeLoc, child)
        typ = return_type(M.grammar, get(child, loc).ind)
        subtree = rand(RuleNode, M.grammar, typ)
        insert!(child, loc, subtree)
    end
    return child
end
####################

#################### expr 5
struct TreePermutation <: MutationMethod
    grammar
	p
end
function mutate(M::TreePermutation, a)
    child = deepcopy(a)
    if rand() < M.p
	    node = sample(child)
        n = length(node.children)
        types = child_types(M.grammar, node)
        for i in 1 : n-1
            c = 1
            for k in i+1 : n
                if types[k] == types[i] &&
                    rand() < 1/(c+=1)

                    node.children[i], node.children[k] =
                        node.children[k], node.children[i]
                end
            end
        end
	end
    return child
end
####################

#################### expr 6
struct DecodedExpression
    node
    n_rules_applied
end
function decode(x, grammar, sym, c_max=1000, c=0)
    node, c = _decode(x, grammar, sym, c_max, c)
    DecodedExpression(node, c)
end
function _decode(x, grammar, typ, c_max, c)
    types = grammar[typ]
    if length(types) > 1
        g = x[mod1(c+=1, length(x))]
        rule = types[mod1(g, length(types))]
    else
        rule = types[1]
    end
    node = RuleNode(rule)
    childtypes = child_types(grammar, node)
    if !isempty(childtypes) && c < c_max
        for ctyp in childtypes
            cnode, c = _decode(x, grammar, ctyp, c_max, c)
            push!(node.children, cnode)
        end
    end
    return (node, c)
end
####################

#################### expr 7
struct IntegerGaussianMutation <: MutationMethod
    σ
end
function mutate(M::IntegerGaussianMutation, child)
    return child + round.(Int, randn(length(child)).*M.σ)
end
####################

#################### expr 8
struct GeneDuplication <: MutationMethod
end
function mutate(M::GeneDuplication, child)
    n = length(child)
    i, j = rand(1:n), rand(1:n)
    interval = min(i,j) : max(i,j)
    return vcat(child, deepcopy(child[interval]))
end
####################

#################### expr 9
struct GenePruning <: MutationMethod
    p
    grammar
    typ
end
function mutate(M::GenePruning, child)
    if rand() < M.p
        c = decode(child, M.grammar, M.typ).n_rules_applied
        if c < length(child)
            child = child[1:c]
        end
    end
    return child
end
####################

#################### expr 10
struct MultiMutate <: MutationMethod
    Ms
end
function mutate(M::MultiMutate, child)
    for m in M.Ms
        child = mutate(m, child)
    end
    return child
end
####################

#################### expr 11
struct ProbabilisticGrammar
    grammar
    ws
end
function probability(probgram, node)
    typ = return_type(probgram.grammar, node)
    i = findfirst(probgram.grammar[typ], node.ind)
    prob = probgram.ws[typ][i] / sum(probgram.ws[typ])
    for (i,c) in enumerate(node.children)
        prob *= probability(probgram, c)
    end
    return prob
end
####################

#################### expr 12
function _update!(probgram, x)
    grammar = probgram.grammar
    typ = return_type(grammar, x)
    i = findfirst(grammar[typ], x.ind)
    probgram.ws[typ][i] += 1
    for c in x.children
        _update!(probgram, c)
    end
    return probgram
end
function update!(probgram, Xs)
    for w in values(probgram.ws)
        fill!(w,0)
    end
    for x in Xs
        _update!(probgram, x)
    end
    return probgram
end
####################

#################### expr 13
struct PPTNode
    ps
    children
end
function PPTNode(grammar;
    w_terminal = 0.6,
    w_nonterm = 1-w_terminal,
    )

    ps = Dict(typ => normalize!([isterminal(grammar, i) ?
                     w_terminal : w_nonterm
                     for i in grammar[typ]], 1)
             for typ in nonterminals(grammar))
    PPTNode(ps, PPTNode[])
end
function get_child(ppt::PPTNode, grammar, i)
    if i > length(ppt.children)
        push!(ppt.children, PPTNode(grammar))
    end
    return ppt.children[i]
end
####################

#################### expr 14
function rand(ppt, grammar, typ)
    rules = grammar[typ]
    rule_index = sample(rules, Weights(ppt.ps[typ]))
    ctypes = child_types(grammar, rule_index)

    arr = Vector{RuleNode}(length(ctypes))
    node = iseval(grammar, rule_index) ?
        RuleNode(rule_index, eval(grammar, rule_index), arr) :
        RuleNode(rule_index, arr)

    for (i,typ) in enumerate(ctypes)
        node.children[i] =
            rand(get_child(ppt, grammar, i), grammar, typ)
    end
    return node
end
####################

#################### expr 15
function probability(ppt, grammar, expr)
    typ = return_type(grammar, expr)
    i = findfirst(grammar[typ], expr.ind)
    prob = ppt.ps[typ][i]
    for (i,c) in enumerate(expr.children)
        prob *= probability(get_child(ppt, grammar, i), grammar, c)
    end
    return prob
end
function p_target(ppt, grammar, x_best, y_best, y_elite, α, ϵ)
    p_best = probability(ppt, grammar, x_best)
    return p_best + (1-p_best)*α*(ϵ - y_elite)/(ϵ - y_best)
end
####################

#################### expr 16
function _update!(ppt, grammar, x, c, α)
    typ = return_type(grammar, x)
    i = findfirst(grammar[typ], x.ind)
    p = ppt.ps[typ]
    p[i] += c*α*(1-p[i])
    psum = sum(p)
    for j in 1 : length(p)
        if j != i
            p[j] *= (1- p[i])/(psum - p[i])
        end
    end
    for (pptchild,xchild) in zip(ppt.children, x.children)
        _update!(pptchild, grammar, xchild, c, α)
    end
    return ppt
end
function update!(ppt, grammar, x_best, y_best, y_elite, α, c, ϵ)
    p_targ = p_target(ppt, grammar, x_best, y_best, y_elite, α, ϵ)
    while probability(ppt, grammar, x_best) < p_targ
        _update!(ppt, grammar, x_best, c, α)
    end
    return ppt
end
####################

#################### expr 17
function mutate!(ppt, grammar, x_best, p_mutation, β;
    sqrtlen = sqrt(length(x_best)),
    )
    typ = return_type(grammar, x_best)
    p = ppt.ps[typ]
    prob = p_mutation/(length(p)*sqrtlen)
    for i in 1 : length(p)
        if rand() < prob
            p[i] += β*(1-p[i])
        end
    end
    normalize!(p, 1)
    for (pptchild,xchild) in zip(ppt.children, x_best.children)
        mutate!(pptchild, grammar, xchild, p_mutation, β,
                sqrtlen=sqrtlen)
    end
    return ppt
end
####################

#################### expr 18
function prune!(ppt, grammar; p_threshold=0.99)
    kmax, pmax = :None, 0.0
    for (k, p) in ppt.ps
        pmax′ = maximum(p)
        if pmax′ > pmax
            kmax, pmax = k, pmax′
        end
    end
    if pmax > p_threshold
        i = indmax(ppt.ps[kmax])
        if isterminal(grammar, i)
            clear!(ppt.children[kmax])
        else
            max_arity_for_rule = maximum(nchildren(grammar, r) for
                                         r in grammar[kmax])
            while length(ppt.children) > max_arity_for_rule
                pop!(ppt.children)
            end
        end
    end
    return ppt
end
####################

#################### mdo 1
function gauss_seidel!(Fs, A; k_max=100, ϵ=1e-4)
	k, converged = 0, false
	while !converged && k ≤ k_max
		k += 1
		A_old = deepcopy(A)
		for F in Fs
			F(A)
		end
		converged = all(isapprox(A[v], A_old[v], rtol=ϵ)
		                for v in keys(A))
	end
	return (A, converged)
end
####################

#################### test-functions 1
function ackley(x, a=20, b=0.2, c=2π)
	d = length(x)
	return -a*exp(-b*sqrt(sum(x.^2)/d)) -
	          exp(sum(cos.(c*xi) for xi in x)/d) + a + e
end
####################

#################### test-functions 2
booth(x) = (x[1]+2x[2]-7)^2 + (2x[1]+x[2]-5)^2
####################

#################### test-functions 3
function branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π))
	return a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
end
####################

#################### test-functions 4
function flower(x; a=1, b=1, c=4)
	return a*norm(x) + b*sin(c*atan2(x[2], x[1]))
end
####################

#################### test-functions 5
function michalewicz(x; m=10)
	return -sum(sin(v)*sin(i*v^2/π)^(2m) for
               (i,v) in enumerate(x))
end
####################

#################### test-functions 6
rosenbrock(x; a=1, b=5) = (a-x[1])^2 + b*(x[2] - x[1]^2)^2
####################

#################### test-functions 7
wheeler(x, a=1.5) = -exp(-(x[1]*x[2] - a)^2 -(x[2]-a)^2)
####################

#################### test-functions 8
function circle(x)
    θ = x[1]
    r = 0.5 + 0.5*(2x[2]/(1+x[2]^2))
    y1 = 1 - r*cos(θ)
    y2 = 1 - r*sin(θ)
    return [y1, y2]
end
####################

#################### math-concepts 1
struct Quadrule
    ws
    xs
end
function quadrule_legendre(m)
    bs = [legendre(i) for i in 1 : m+1]
    xs = roots(bs[end])
    A = [bs[k](xs[i]) for k in 1 : m, i in 1 : m]
    b = zeros(m)
    b[1] = 2
    ws = A\b
    return Quadrule(ws, xs)
end
####################

#################### math-concepts 2
quadint(f, quadrule) =
    sum(w*f(x) for (w,x) in zip(quadrule.ws, quadrule.xs))
function quadint(f, quadrule, a, b)
    α = (b-a)/2
    β = (a+b)/2
    g = x -> α*f(α*x+β)
    return quadint(g, quadrule)
end
####################

