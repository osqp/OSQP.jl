using OSQP, Base.Test


function setup()
        options = Dict(:verbose => false,
                       :eps_abs => 1e-06,
                       :eps_rel => 1e-06,
                       :auto_rho => false,
                       :max_iter => 5000)
	return options
end

tol = 1e-3


@testset "feasibility" begin

	@testset "feasibility_problem" begin
		n = 30
		m = 30
		P = spzeros(n, n) 
		q = zeros(n)
		A = sprandn(m, n, 0.8) 
		u = randn(m) 
		l = u 
		options = setup()

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
		results = OSQP.solve!(model)	

		@test isapprox(norm(A * results.x - u), 0., atol=tol) 


	end

end
