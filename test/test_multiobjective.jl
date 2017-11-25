let
	@test !dominates([1.0, 2.0], [1.0, 2.0])
	@test  dominates([1.0, 2.0], [1.0, 2.5])
	@test  dominates([1.0, 2.0], [1.5, 2.0])
	@test  dominates([-1.0, -2.0], [-1.0, -1.5])

	xs = [
		[1.0],[2.0],[3.0],[4.0]
	]
	ys = [
		[0.0,0.0],[2.0,-1.0],[-1.0,2.0],[2.0,2.0]
	]
	pareto_xs, pareto_ys = naive_pareto(xs, ys)
	@test pareto_xs == [[1.0],[2.0],[3.0]]
	@test pareto_ys == [[0.0,0.0],[2.0,-1.0],[-1.0,2.0]]

	xs = Vector{Float64}[]
	ys = Vector{Float64}[]
	for i in 1 : 100
	    push!(xs, rand(2))
	    push!(ys, rand(2))
	end
	naive_pareto(xs, ys)

	xs = weight_pareto(x->-cos(x), x->-sin(x), 5)
	@test isapprox(xs[5], atan2(0,4), atol=1e-6)
	@test isapprox(xs[4], atan2(1,3), atol=1e-6)
	@test isapprox(xs[3], atan2(2,2), atol=1e-6)
	@test isapprox(xs[2], atan2(3,1), atol=1e-6)
	@test isapprox(xs[1], atan2(4,0), atol=1e-6)

	population = [randn(1) for i in 1 : 20]
	vector_evaluated_genetic_algorithm(x->[-cos(x[1]), -sin(x[1])], population, 10, TruncationSelection(5), SinglePointCrossover(), GaussianMutation(0.01))

	get_non_domination_levels(ys)

	srand(0)
	xs = [
		[1.0],[2.0],[3.0],[4.0]
	]
	ys = [
		[0.0,0.0],[2.0,-1.0],[-2.0,2.0],[2.0,2.0]
	]
	(xs, ys) = discard_closest_pair!(xs, ys)
	@test (xs == [[1.0],[3.0],[4.0]] && ys == [[0.0, 0.0],[-2.0,2.0],[2.0,2.0]]) ||
	      (xs == [[2.0],[3.0],[4.0]] && ys == [[2.0,-1.0],[-2.0,2.0],[2.0,2.0]])

	filter_xs = Vector{Float64}[]
	filter_ys = Vector{Float64}[]
	for i in 1 : 10
	    push!(filter_xs, rand(2))
	    push!(filter_ys, rand(2))
	end
	update_pareto_filter!(filter_xs, filter_ys, xs, ys)
end