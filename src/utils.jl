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
