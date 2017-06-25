# Wrapper for the low level functions defined in https://github.com/oxfordcontrol/osqp/blob/master/include/osqp.h


# macro to call a OSQP C function
# macro osqp_ccall(func, args...)
#     f = "osqp_$(func)"
#     args = map(esc,args)
#     is_unix() && return quote
#         ccall(($f,OSQP.osqp), $(args...))
#     end
#     is_windows() && VERSION < v"0.6-" && return quote
#         ccall(($f,OSQP.osqp), stdcall, $(args...))
#     end
#     is_windows() && VERSION >= v"0.6-" && return quote
#         ccall(($f,OSQP.osqp), $(esc(:stdcall)), $(args...))
#     end
# end
#

type Model 
	"""OSQP workspace pointer"""
	workspace::Ptr{OSQP.Workspace}

	"""
	    setup(P, a, A, l, u, settings)
	    
	Perform OSQP solver setup using the inputs `P`, `q`, `A`, `l`, `u`
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

function setup!(model::OSQP.Model,
		P::SparseMatrixCSC=nothing, 
		q::Vector{Float64}=nothing, 
		A::SparseMatrixCSC=nothing, 
		l::Vector{Float64}=nothing,
		u::Vector{Float64}=nothing;
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

	# Do not use this anymore. We instead copy the solution
	# x = unsafe_wrap(Array, solution.x, data.n)
	# y = unsafe_wrap(Array, solution.y, data.m)

	# Allocate solution vectors and copy solution
	x = Array{Float64}(data.n)
	y = Array{Float64}(data.m)
	unsafe_copy!(pointer(x), solution.x, data.n)
	unsafe_copy!(pointer(y), solution.y, data.m)
	
	# Recover Cinfo structure
        cinfo = unsafe_load(workspace.info)	

	# Construct C structure
	info = OSQP.Info(cinfo)

	# Return results
	return Results(x, y, info)


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
			if !(key in updatable_data)
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
			if !(key in updatable_settings)
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
	alpha = get(data, :alpha, nothing)
	delta = get(data, :delta, nothing)
	polish = get(data, :polish, nothing)
	pol_refine_iter = get(data, :pol_refine_iter, nothing)
	verbose = get(data, :verbose, nothing)
	early_terminate = get(data, :early_terminate, nothing)
	early_terminate_interval = get(data, :early_terminate_interval, nothing)
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

	if pol_refine_iter != nothing
		exitflag = ccall((:osqp_update_pol_refine_iter, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, pol_refine_iter) 
		if exitflag != 0 error("Error updating pol_refine_iter") end
	end

	if verbose != nothing
		exitflag = ccall((:osqp_update_verbose, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, verbose) 
		if exitflag != 0 error("Error updating verbose") end
	end

	if early_terminate != nothing
		exitflag = ccall((:osqp_update_early_terminate, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, early_terminate) 
		if exitflag != 0 error("Error updating early_terminate") end
	end

	if early_terminate_interval != nothing
		exitflag = ccall((:osqp_update_early_terminate_interval, OSQP.osqp), Clong, (Ptr{OSQP.Workspace}, Clong), model.workspace, early_terminate_interval) 
		if exitflag != 0 error("Error updating early_terminate_interval") end
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
