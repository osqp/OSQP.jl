# Wrapper for the low level functions defined in https://github.com/oxfordcontrol/osqp/blob/master/include/osqp.h

type Model
	workspace::Ptr{OSQP.Workspace}

	"""
	    Module()

	Initialize OSQP module
	"""
	function Model()
			# TODO: Change me to more elegant way
			# a = Array{Ptr{OSQP.Workspace}}(1)[1]
			a = C_NULL

			# Create new model
			model = new(a)

			# Add finalizer
			finalizer(model, OSQP.clean!)

			return model

	end


end

"""
    setup!(module, P, a, A, l, u, settings)

Perform OSQP solver setup of module `module`, using the inputs `P`, `q`, `A`, `l`, `u`
"""
function setup!(model::OSQP.Model;
		P::Union{SparseMatrixCSC, Void}=nothing,
		q::Union{Vector{Float64}, Void}=nothing,
		A::Union{SparseMatrixCSC, Void}=nothing,
		l::Union{Vector{Float64}, Void}=nothing,
		u::Union{Vector{Float64}, Void}=nothing,
		settings...)

	# Check problem dimensions
	if P == nothing
		if q != nothing
			n = length(q)
		elseif A != nothing
			n = size(A, 2)
		else
			error("The problem does not have any variables!")
		end

	else
		n = size(P, 1)
	end

	if A == nothing
		m = 0
	else
		m = size(A, 1)
	end


	# Check if parameters are nothing
	if ((A == nothing) & ( (l != nothing) | (u != nothing))) |
		((A != nothing) & ((l == nothing) | (u == nothing)))
		error("A must be supplied together with l and u")
	end

	if (A != nothing) & (l == nothing)
		l = - Inf * ones(m)
	end
	if (A != nothing) & (u == nothing)
		u = Inf * ones(m)
	end

	if P == nothing
		P = sparse([], [], [], n, n)
	end
	if q == nothing
		q = zeros(n)
	end
	if A == nothing
		A = sparse([], [], [], m, n)
		l = zeros(m)
		u = zeros(m)
	end


	# Check if dimensions are correct
	if length(q) != n
		error("Incorrect dimension of q")
	end
	if length(l) != m
		error("Incorrect dimensions of l")
	end
	if length(u) != m
		error("Incorrect dimensions of u")
	end


	# Check or sparsify matrices
	if !issparse(P)
		warn("P is not sparse. Sparsifying it now (it might take a while)")
		P = sparse(P)
	end
	if !issparse(A)
		warn("A is not sparse. Sparsifying it now (it might take a while)")
		A = sparse(A)
	end

	# Convert lower and upper bounds from Julia infinity to OSQP infinity
	u = min.(u, OSQP_INFTY)
	l = max.(l, -OSQP_INFTY)

	# Create OSQP data
	data = OSQP.Data(n, m, P, q, A, l, u)

	# Create OSQP settings
	stgs = OSQP.Settings(settings)

	# Perform setup
	model.workspace = ccall((:osqp_setup, OSQP.osqp),
				Ptr{OSQP.Workspace}, (Ptr{OSQP.Data},
	                        Ptr{OSQP.Settings}), &data, &stgs)

	if model.workspace == C_NULL
		error("Error in OSQP setup")
	end

end


function solve!(model::OSQP.Model)

	# Solve problem
	exitflag = ccall((:osqp_solve, OSQP.osqp), Clong,
			 (Ptr{OSQP.Workspace}, ),
			 model.workspace)

	if exitflag != 0
		error("Error in OSQP solution!")
	end

	# Recover solution
	workspace = unsafe_load(model.workspace)
	solution = unsafe_load(workspace.solution)
	data = unsafe_load(workspace.data)

	# Recover Cinfo structure
        cinfo = unsafe_load(workspace.info)

	# Construct C structure
	info = OSQP.Info(cinfo)

	# Do not use this anymore. We instead copy the solution
	# x = unsafe_wrap(Array, solution.x, data.n)
	# y = unsafe_wrap(Array, solution.y, data.m)

	# Allocate solution vectors and copy solution
	x = Array{Float64}(data.n)
	y = Array{Float64}(data.m)

        if info.status in SOLUTION_PRESENT 
		# If solution exists, copy it
		unsafe_copy!(pointer(x), solution.x, data.n)
		unsafe_copy!(pointer(y), solution.y, data.m)

		# Return results
		return Results(x, y, info)
	else
		# else fill with NaN and return certificates of infeasibility
		x *= NaN
		y *= NaN
		if info.status == :Primal_infeasible || info.status == :Primal_infeasible_inaccurate
			prim_inf_cert = Array{Float64}(data.m) 
			unsafe_copy!(pointer(prim_inf_cert), workspace.delta_y, data.m)
			# Return results
		        return Results(x, y, info, prim_inf_cert, nothing)	
		elseif info.status == :Dual_infeasible || info.status == :Dual_infeasible_inaccurate
			dual_inf_cert = Array{Float64}(data.n) 
			unsafe_copy!(pointer(dual_inf_cert), workspace.delta_x, data.n)
			# Return results
		        return Results(x, y, info, nothing, dual_inf_cert)	
		end
	end
end


function version()
    return unsafe_string(ccall((:osqp_version, OSQP.osqp), Cstring, ()))
end

function clean!(model::OSQP.Model)
	exitflag = ccall((:osqp_cleanup, OSQP.osqp), Clong,
			 (Ptr{OSQP.Workspace},), model.workspace)
	if exitflag != 0
		error("Error in OSQP cleanup")
	end
end


function update!(model::OSQP.Model; kwargs...)
	if isempty(kwargs)
		return
	else
		data = Dict{Symbol, Any}()
		for (key, value) in kwargs
			if !(key in UPDATABLE_DATA)
				error("$(key) field cannot be updated or is not recognized")
			else
				data[key] = value
			end
		end
	end

	# Get arguments
	q = get(data, :q, nothing)
	l = get(data, :l, nothing)
	u = get(data, :u, nothing)
	Px = get(data, :Px, nothing)
	Px_idx = get(data, :Px_idx, C_NULL)
	Ax = get(data, :Ax, nothing)
	Ax_idx = get(data, :Ax_idx, C_NULL)

	# Get problem dimensions
	(n, m) = OSQP.dimensions(model)

	# Update linear cost
	if q != nothing
		if length(q) != n
			error("q must have length n = $(n)")
		end
		exitflag = ccall((:osqp_update_lin_cost, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}), model.workspace, pointer(q))
		if exitflag != 0 error("Error updating q") end
	end


	# Update lower bound
	if l != nothing
		if length(l) != m
			error("l must have length m = $(m)")
		end


		# Convert values to OSQP_INFTY
		l = max.(l, -OSQP_INFTY)

		if u == nothing
			exitflag = ccall((:osqp_update_lower_bound, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}), model.workspace, pointer(l))
			if exitflag != 0 error("Error updating l") end
		end
	end


	# Update upper bound
	if u != nothing
		if length(u) != m
			error("u must have length m = $(m)")
		end


		# Convert values to OSQP_INFTY
		u = min.(u, OSQP_INFTY)

		if l == nothing
			exitflag = ccall((:osqp_update_upper_bound, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}), model.workspace, pointer(u))
			if exitflag != 0 error("Error updating u") end
		end
	end


	# Update bounds
	if (l != nothing) & (u != nothing)
		exitflag = ccall((:osqp_update_bounds, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}, Ptr{Cdouble}), model.workspace, pointer(l), pointer(u))
		if exitflag != 0 error("Error updating bounds l and u") end
	end


	# Update matrix P
	if Px != nothing
		if (Px_idx != C_NULL)
			if (length(Px_idx) != length(Px))
				error("Px and Px_idx must have same length")
			end
			Px_idx = pointer(Px_idx)  # Get pointer to pass to the C function
		end
		if Ax == nothing
			exitflag = ccall((:osqp_update_P, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}, Ptr{Clong}, Clong),
					 model.workspace, pointer(Px), Px_idx, length(Px))
			if exitflag != 0 error("Error updating P") end
		end
	end

	# Update matrix A
	if Ax != nothing
		if (Ax_idx != C_NULL)
			if (length(Ax_idx) != length(Ax))
				error("Ax and Ax_idx must have same length")
			end
			Ax_idx = pointer(Ax_idx)  # Get pointer to pass to the C function
		end
		if Px == nothing
			exitflag = ccall((:osqp_update_A, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}, Ptr{Clong}, Clong),
					 model.workspace, pointer(Ax), Ax_idx, length(Ax))
			if exitflag != 0 error("Error updating A") end
		end
	end

	# Update both matrices P and A
	if (Px != nothing) & (Ax != nothing)
		exitflag = ccall((:osqp_update_P_A, OSQP.osqp), Clong,
				 (Ptr{OSQP.Workspace}, Ptr{Cdouble}, Ptr{Clong}, Clong,
				  Ptr{Cdouble}, Ptr{Clong}, Clong),
				 model.workspace, pointer(Px), Px_idx, length(Px), pointer(Ax), Ax_idx, length(Ax))
		if exitflag != 0 error("Error updating P and A") end
	end

end



function update_settings!(model::OSQP.Model; kwargs...)

	if isempty(kwargs)
		return
	else
		data = Dict{Symbol, Any}()
		for (key, value) in kwargs
			if !(key in UPDATABLE_SETTINGS)
				error("$(key) cannot be updated or is not recognized")
			else
				data[key] = value
			end
		end
	end

	# Get arguments
	max_iter = get(data, :max_iter, nothing)
	eps_abs = get(data, :eps_abs, nothing)
	eps_rel = get(data, :eps_rel, nothing)
	eps_prim_inf = get(data, :eps_prim_inf, nothing)
	eps_dual_inf = get(data, :eps_dual_inf, nothing)
	rho = get(data, :rho, nothing)
	alpha = get(data, :alpha, nothing)
	delta = get(data, :delta, nothing)
	polish = get(data, :polish, nothing)
	polish_refine_iter = get(data, :polish_refine_iter, nothing)
	verbose = get(data, :verbose, nothing)
	scaled_termination = get(data, :early_terminate, nothing)
	check_termination = get(data, :check_termination, nothing)
	warm_start = get(data, :warm_start, nothing)

	# Update individual settings
	if max_iter != nothing
		exitflag = ccall((:osqp_update_max_iter, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, max_iter)
		if exitflag != 0 error("Error updating max_iter") end
	end

	if eps_abs != nothing
		exitflag = ccall((:osqp_update_eps_abs, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, eps_abs)
		if exitflag != 0 error("Error updating eps_abs") end
	end

	if eps_rel != nothing
		exitflag = ccall((:osqp_update_eps_rel, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, eps_rel)
		if exitflag != 0 error("Error updating eps_rel") end
	end


	if eps_prim_inf != nothing
		exitflag = ccall((:osqp_update_eps_prim_inf, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, eps_prim_inf)
		if exitflag != 0 error("Error updating eps_prim_inf") end
	end

	if eps_dual_inf != nothing
		exitflag = ccall((:osqp_update_eps_dual_inf, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, eps_dual_inf)
		if exitflag != 0 error("Error updating eps_dual_inf") end
	end

	if rho != nothing
		exitflag = ccall((:osqp_update_rho, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, rho)
		if exitflag != 0 error("Error updating rho") end
	end

	if alpha != nothing
		exitflag = ccall((:osqp_update_alpha, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, alpha)
		if exitflag != 0 error("Error updating alpha") end
	end

	if delta != nothing
		exitflag = ccall((:osqp_update_delta, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Cdouble), model.workspace, delta)
		if exitflag != 0 error("Error updating delta") end
	end

	if polish != nothing
		exitflag = ccall((:osqp_update_polish, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, polish)
		if exitflag != 0 error("Error updating polish") end
	end

	if polish_refine_iter != nothing
		exitflag = ccall((:osqp_update_polish_refine_iter, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, polish_refine_iter)
		if exitflag != 0 error("Error updating polish_refine_iter") end
	end

	if verbose != nothing
		exitflag = ccall((:osqp_update_verbose, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, verbose)
		if exitflag != 0 error("Error updating verbose") end
	end

	if scaled_termination != nothing
		exitflag = ccall((:osqp_update_scaled_termination, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, scaled_termination)
		if exitflag != 0 error("Error updating scaled_termination") end
	end

	if check_termination != nothing
		exitflag = ccall((:osqp_update_check_termination, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, check_termination)
		if exitflag != 0 error("Error updating check_termination") end
	end

	if warm_start != nothing
		exitflag = ccall((:osqp_update_warm_start, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, warm_start)
		if exitflag != 0 error("Error updating warm_start") end
	end

	return nothing
end



function warm_start!(model::OSQP.Model; x::Vector{Float64}=nothing, y::Vector{Float64}=nothing)
	# Get problem dimensions
	(n, m) = OSQP.dimensions(model)

	if x != nothing
		if length(x) != n
			error("Wrong dimension for variable x")
		end

		if y == nothing
			exitflag = ccall((:osqp_warm_start_x, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}), model.workspace, x)
			if exitflag != 0 error("Error in warm starting x") end
		end
	end


	if y != nothing
		if length(y) != m
			error("Wrong dimension for variable y")
		end

		if x == nothing
			exitflag = ccall((:osqp_warm_start_y, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}), model.workspace, y)
			if exitflag != 0 error("Error in warm starting y") end
		end
	end

	if (x != nothing) & (y != nothing)
		exitflag = ccall((:osqp_warm_start, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Ptr{Cdouble}, Ptr{Cdouble}), model.workspace, x, y)
		if exitflag != 0 error("Error in warm starting x and y") end
	end

end



# Auxiliary low-level functions
"""
    dimensions(model::OSQP.Model)

Obtain problem dimensions from OSQP model
"""
function dimensions(model::OSQP.Model)

	workspace = unsafe_load(model.workspace)
	if workspace == C_NULL
		error("Workspace has not been setup yet")
	end
	data = unsafe_load(workspace.data)
	return data.n, data.m
end




function linsys_solver_str_to_int!(settings_dict::Dict{Symbol, Any})
         # linsys_str = pop!(settings_dict, :linsys_solver)
         linsys_str = get(settings_dict, :linsys_solver, nothing)	

	 if linsys_str != nothing
		 # Check type
		 if !isa(linsys_str, String)
			error("linsys_solver is required to be a string")
		 end

		 # Convert to lower case
		 linsys_str = lowercase(linsys_str)

		 if linsys_str == "suitesparse ldl"
			settings_dict[:linsys_solver] = SUITESPARSE_LDL_SOLVER
		elseif linsys_str == "mkl pardiso"
			settings_dict[:linsys_solver] = MKL_PARDISO_SOLVER
		elseif linsys_str == ""
			settings_dict[:linsys_solver] = SUITESPARSE_LDL_SOLVER
		else
			warn("Linear system solver not recognized. Using default SuiteSparse LDL")
			settings_dict[:linsys_solver] = SUITESPARSE_LDL_SOLVER

		end	
	end
end
