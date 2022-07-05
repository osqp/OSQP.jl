const OSQP_UNKNOWN_SOLVER = 0
const OSQP_DIRECT_SOLVER = 1
const OSQP_INDIRECT_SOLVER = 2

# Define OSQP infinity constants
const OSQP_INFTY = 1e30

# OSQP return values
# https://github.com/oxfordcontrol/osqp/blob/master/include/constants.h
const status_map = Dict{Int,Symbol}(
    1 => :Solved,
    2 => :Solved_inaccurate,
    3 => :Primal_infeasible,
    4 => :Primal_infeasible_inaccurate,
    5 => :Dual_infeasible,
    6 => :Dual_infeasible_inaccurate,
    7 => :Max_iter_reached,
    8 => :Time_limit_reached,
    9 => :Non_convex,
    10 => :Interrupted,
    11 => :Unsolved,
)

const SOLUTION_PRESENT = [:Solved_inaccurate, :Solved, :Max_iter_reached]

# UPDATABLE_DATA
const UPDATABLE_DATA = [:q, :l, :u, :Px, :Px_idx, :Ax, :Ax_idx]
