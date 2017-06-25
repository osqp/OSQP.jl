# Utilities
# N.B. Cannot call constants directly. Constants in headers are not shared in dynamically linked library
# function constant(constant_name)
#         constant = ccall((:osqp_constant, OSQP.osqp), Clong, (Ptr{Clong}))
#         return constant
# end



# Define OSQP infinity constants
OSQP_INFTY = 1e20

# OSQP return values
# https://github.com/oxfordcontrol/osqp/blob/master/include/constants.h
const status_map = Dict{Int, Symbol}(
    1 => :Solved,
    -2 => :Max_Iter_Reached,
    -3 => :Primal_Infeasible,
    -4 => :Dual_Infeasible,
    -5 => :Interrupted,
    -10 => :Unsolved
)


# updatable_data
updatable_data = [:q, :l, :u, :Px, :Px_idx, :Ax, :Ax_idx]

# updatable_settings
updatable_settings = [:max_iter, :eps_aps, :eps_rel, :eps_prim_inf, :eps_dual_inf,
		      :alpha, :delta, :polish, :pol_refine_iter, :verbose, :early_terminate,
		      :early_terminate_interval, :warm_start]


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
