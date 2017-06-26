using OSQP, Base.Test


function setup()
        options = Dict(:verbose => false,
                       :eps_abs => 1e-08,
                       :eps_rel => 1e-08,
                       :polish => false,
		       :auto_rho => false,
		       :early_terminate_interval => 1)
	return options
end

tol = 1e-5


@testset "warm_start" begin

	@testset "warm_start_problem" begin
		srand(1)

		n = 100 
		m = 200
		P = sprandn(n, n, 0.9)
		P = P' * P
		q = randn(n) 
		A = sprandn(m, n, 0.9)
		u = rand(m) * 2
		l = -rand(m) * 2



		options = setup()

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
		results = OSQP.solve!(model)	

		# Store optimal values
		x_opt = results.x
		y_opt = results.y
		tot_iter = results.info.iter


		# Warm start with zeros to check if number of iterations is the same
		OSQP.warm_start!(model, x=zeros(n), y=zeros(m))
		results = OSQP.solve!(model)
		@test results.info.iter == tot_iter
		
		# Warm start with optimal values and check that total number of iter < 10
		OSQP.warm_start!(model, x=x_opt, y=y_opt)
		results = OSQP.solve!(model)
		@test results.info.iter <= 10 


	end

end
