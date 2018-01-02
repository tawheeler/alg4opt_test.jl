let
	# find (x₁, x₂) such that x₂ = 2*x₁ and x₁ = x₂ - 1
	# has solution [1,2]

	sol1 = Dict(:x₁=>0.0, :x₂=>0.0)
	sol2 = Dict(:x₁=>2.0, :x₂=>4.0)

	function F1(A)
		x₁ = A[:x₁]
		A[:x₂] = 2x₁
		return A
	end
	function F2(A)
		x₂ = A[:x₂]
		A[:x₁] = sqrt(x₂)
		return A
	end
	Fs = [F1,F2]

	# these should converge right away
	@test gauss_seidel!(Fs, sol1)[2]
	@test gauss_seidel!(Fs, sol2)[2]

	A, converged = gauss_seidel!(Fs, Dict(:x₁=>1.6, :x₂=>NaN), ϵ=1e-5)
	@test converged
	@test isapprox(A[:x₁], sol2[:x₁], atol=1e-4)
	@test isapprox(A[:x₂], sol2[:x₂], atol=1e-4)
end
let
	# this one does not converge
	function F1(A)
		x₁ = A[:x₁]
		A[:x₂] = 2x₁
		return A
	end
	function F2(A)
		x₂ = A[:x₂]
		A[:x₁] = x₂^2
		return A
	end
	Fs = [F1,F2]

	A, converged = gauss_seidel!(Fs, Dict(:x₁=>1.6, :x₂=>NaN), k_max=5)
	@show (A, converged)
	@test !converged
end