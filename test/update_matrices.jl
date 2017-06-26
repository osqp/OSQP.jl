using OSQP, Base.Test


function setup()
        options = Dict(:verbose => false,
                       :eps_abs => 1e-08,
                       :eps_rel => 1e-08,
                       :polish => false,
		       :auto_rho => false,
		       :early_terminate_interval => 1)
	
	srand(1)

	n = 5
	m = 8
	p = 0.7
	Pt = sprandn(n, n, p)
	Pt_new = copy(Pt)
	P = Pt * Pt' + speye(n)
	(Pti, Ptj) = findn(Pt)
	Ptx = copy(Pt.nzval)
	
	Pt_newx = Ptx + 0.1 * randn(length(Ptx))
	# Pt_new = sparse(Pi, Pj, Pt_newx)
	P_new = Pt_new * Pt_new' + speye(n)
	q = randn(n)
	A = sprandn(m, n, p)


	(Ai, Aj) = findn(A)
	Ax = copy(A.nzval)
	A_newx = Ax + randn(length(Ax))
	A_new = sparse(Ai, Aj, A_newx)

	# A_new = copy(A)
	# A_new.nzval += randn(length(A_new.nzval))
	l = zeros(m)
	u = 30 + randn(m)

	
	problem = Dict()
	problem[:P] = P
	problem[:P_new] = P_new
	problem[:q] = q
	problem[:A] = A
	problem[:A_new] = A_new
	problem[:l] = l
	problem[:u] = u
	problem[:m] = m
	problem[:n] = n
	return problem, options
end

tol = 1e-5

@testset "update_matrices" begin

	@testset "solve" begin
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
		results = OSQP.solve!(model)	

		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")


		x_test = [-0.324865, 0.598681, -0.066646, -0.00653471, -0.0736556]
		y_test = [-2.72542e-10, -1.927e-12, -0.547374, -2.25907e-10, -0.645086, -4.09874e-9, -1.4083e-11, -1.26234e-11]
		obj_test = -0.6736977821770016
		
		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 


	end


        @testset "update_P" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrix
		Pnew_triu = triu(P_new)
		# Note the shift in indexing! (Cherent with C)
		Pnew_triu_idx = collect(range(1, length(Pnew_triu.nzval))) - 1
		OSQP.update!(model, Px=Pnew_triu.nzval, Px_idx=Pnew_triu_idx)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")

		x_test = [-0.324865, 0.598681, -0.066646, -0.00653471, -0.0736556]
		y_test = [-2.72542e-10, -1.927e-12, -0.547374, -2.25907e-10, -0.645086, -4.09874e-9, -1.4083e-11, -1.26234e-11]
		obj_test = -0.6736977821770016

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_P_allind" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrix
		Pnew_triu = triu(P_new)
		OSQP.update!(model, Px=Pnew_triu.nzval)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")

		x_test = [-0.324865, 0.598681, -0.066646, -0.00653471, -0.0736556]
		y_test = [-2.72542e-10, -1.927e-12, -0.547374, -2.25907e-10, -0.645086, -4.09874e-9, -1.4083e-11, -1.26234e-11]
		obj_test = -0.6736977821770016

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_A" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrix
		# Note the shift in indexing! (Cherent with C)
		A_new_idx = collect(range(1, length(A_new.nzval))) - 1
		OSQP.update!(model, Ax=A_new.nzval, Ax_idx=A_new_idx)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703


		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

	@testset "update_A_allind" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrix
		OSQP.update!(model, Ax=A_new.nzval)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703


		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_P_A_indP_indA" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrices P and A
		P_new_triu = triu(P_new)
		# Note the shift in indexing! (Cherent with C)
		P_new_triu_idx = collect(range(1, length(P_new_triu.nzval))) - 1
		A_new_idx = collect(range(1, length(A_new.nzval))) - 1

		OSQP.update!(model, Px=P_new_triu.nzval, Px_idx=P_new_triu_idx,
			     Ax=A_new.nzval, Ax_idx=A_new_idx)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_P_A_indP" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrices P and A
		P_new_triu = triu(P_new)
		# Note the shift in indexing! (Cherent with C)
		P_new_triu_idx = collect(range(1, length(P_new_triu.nzval))) - 1

		OSQP.update!(model, Px=P_new_triu.nzval, Px_idx=P_new_triu_idx,
			     Ax=A_new.nzval)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_P_A_indA" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrices P and A
		P_new_triu = triu(P_new)
		# Note the shift in indexing! (Cherent with C)
		A_new_idx = collect(range(1, length(A_new.nzval))) - 1

		OSQP.update!(model, Px=P_new_triu.nzval,
			     Ax=A_new.nzval, Ax_idx=A_new_idx)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

        @testset "update_P_A_allind" begin
		
		problem, options = setup()
		
		(n, m, P, q, A, l, u) = (problem[:n], problem[:m], 
					 problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
		P_new = problem[:P_new]
		A_new = problem[:A_new]

		model = OSQP.Model()
		OSQP.setup!(model, P, q, A, l, u; options...)
	
		# Update matrices P and A
		P_new_triu = triu(P_new)

		OSQP.update!(model, Px=P_new_triu.nzval, Ax=A_new.nzval)
		results = OSQP.solve!(model)
		
		# # Solve with Gurobi
		# using Gurobi
		# env = Gurobi.Env()
		# setparam!(env, "OutputFlag", 1)
		# model = Gurobi.Model(env, "model")
		# add_cvars!(model, q); update_model!(model)
		# add_constrs!(model, A_new, repmat(['<'], m), u); update_model!(model)
		# add_constrs!(model, A_new, repmat(['>'], m), l); update_model!(model)
		# add_qpterms!(model, P_new); update_model!(model)
		# optimize(model)
		# x_test = get_solution(model)
		# obj_test= get_objval(model)
		# y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
		# y_test = y_test[m+1:end] + y_test[1:m]
                #
		# println("x_test = $(x_test)")
		# println("y_test = $(y_test)")
		# println("obj_test = $(obj_test)")
                #
		x_test = [-0.0318085, 0.067069, -0.0242966, -0.0593736, -0.0321274]
		y_test = [-2.16989e-10, -1.01848, -0.516013, -0.0227263, -3.19721e-10, -1.04482, -3.4596e-10, -4.51608e-10]
		obj_test = -0.02187014865840703

		@test isapprox(results.x, x_test, atol=tol)
		@test isapprox(results.y, y_test, atol=tol)
		@test isapprox(results.info.obj_val, obj_test, atol=tol)
		@test results.info.status == :Solved 
	end

end
