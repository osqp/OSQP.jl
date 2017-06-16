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
			a = Array{Ptr{OSQP.Workspace}}(1)[1]
			
			# Create new model 
			model = new(a)

			# Add finalizer
			finalizer(model, OSQP.clean)

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
	settings = OSQP.Settings()

	# Perform setup
	model.workspace = ccall((:osqp_setup, OSQP.osqp), 
				Ptr{OSQP.Workspace}, (Ptr{OSQP.Data}, 
	                        Ptr{OSQP.Settings}), &data, &settings)


end


function solve(model::OSQP.Model)
	
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
	x = unsafe_wrap(Array, solution.x, data.n)
	y = unsafe_wrap(Array, solution.y, data.m)

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

function clean(model::OSQP.Model)
	exitflag = ccall((:osqp_cleanup, OSQP.osqp), Clong, 
			 (Ptr{OSQP.Workspace},), model.workspace)
	if exitflag != 0
		error("Error in OSQP cleanup")
	end
end


# TODO: Use get/set default settings from C to set default settings in setup
