# Define all the functions in the OSQP C library in a single place
# Tint represents the integer type of the underlying C library
# Tfloat represents the float type of the underlying library
@cprototype osqp_capabilities()::Tint

@cprototype osqp_version()::Cstring

@cprototype osqp_error_message(error_flag::Tint)::Cstring

@cprototype osqp_get_dimensions(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, m::Ptr{Tint}, n::Ptr{Tint})::Nothing

@cprototype osqp_set_default_settings(settings::Ref{OSQP.OSQPSettings})::Nothing

@cprototype osqp_setup(solverp::Ptr{Ptr{OSQP.OSQPSolver{Tfloat,Tint}}}, P::Ptr{OSQP.Ccsc{Tfloat,Tint}}, q::Ptr{Tfloat}, A::Ptr{OSQP.Ccsc{Tfloat,Tint}}, l::Ptr{Tfloat}, u::Ptr{Tfloat}, m::Tint, n::Tint, settings::Ptr{OSQP.OSQPSettings{Tfloat,Tint}})::Tint

@cprototype osqp_solve(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}})::Tint

@cprototype osqp_cleanup(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}})::Tint

@cprototype osqp_warm_start(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, x::Ptr{Tfloat}, y::Ptr{Tfloat})::Tint

@cprototype osqp_cold_start(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}})::Nothing

@cprototype osqp_update_data_vec(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, q_new::Ptr{Tfloat}, l_new::Ptr{Tfloat}, u_new::Ptr{Tfloat})::Tint

@cprototype osqp_update_data_mat(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, Px_new::Ptr{Tfloat}, Px_new_idx::Ptr{Tint}, P_new_n::Tint, Ax_new::Ptr{Tfloat}, Ax_new_idx::Ptr{Tint}, A_new_n::Tint)::Tint

@cprototype osqp_update_settings(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, settings::Ptr{OSQP.OSQPSettings{Tfloat,Tint}})::Tint

@cprototype osqp_update_rho(solver::Ptr{OSQP.OSQPSolver{Tfloat,Tint}}, rho_new::Tfloat)::Tint
