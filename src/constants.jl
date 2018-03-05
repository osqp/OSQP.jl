const SUITESPARSE_LDL_SOLVER = 0
const MKL_PARDISO_SOLVER = 1

# Define OSQP infinity constants
const OSQP_INFTY = 1e20

# OSQP return values
# https://github.com/oxfordcontrol/osqp/blob/master/include/constants.h
const status_map = Dict{Int,Symbol}(4 => :Dual_infeasible_inaccurate,
    3 => :Primal_infeasible_inaccurate,
    2 => :Solved_inaccurate,
    1 => :Solved,
    -2 => :Max_iter_reached,
    -3 => :Primal_infeasible,
    -4 => :Dual_infeasible,
    -5 => :Interrupted,
    -6 => :Time_limit_reached,
    -10 => :Unsolved)

const SOLUTION_PRESENT = [:Solved_inaccurate, :Solved, :Max_iter_reached]

# UPDATABLE_DATA
const UPDATABLE_DATA = [:q, :l, :u, :Px, :Px_idx, :Ax, :Ax_idx]

# UPDATABLE_SETTINGS
const UPDATABLE_SETTINGS = [:max_iter, :eps_aps, :eps_rel, :eps_prim_inf, :eps_dual_inf, :time_limit,
              :rho, :alpha, :delta, :polish, :polish_refine_iter, :verbose,
              :check_termination,:warm_start]


