let
	@test !dominates([1.0, 2.0], [1.0, 2.0])
	@test  dominates([1.0, 2.0], [1.0, 2.5])
	@test  dominates([1.0, 2.0], [1.5, 2.0])
	@test  dominates([-1.0, -2.0], [-1.0, -1.5])

	xs = Vector{Float64}[]
	ys = Vector{Float64}[]
	for i in 1 : 100
	    push!(xs, rand(2))
	    push!(ys, rand(2))
	end
	naive_pareto(xs, ys)

	get_non_domination_levels(ys)

	filter_xs = Vector{Float64}[]
	filter_ys = Vector{Float64}[]
	for i in 1 : 10
	    push!(filter_xs, rand(2))
	    push!(filter_ys, rand(2))
	end
	update_pareto_filter!(filter_xs, filter_ys, xs, ys)
end