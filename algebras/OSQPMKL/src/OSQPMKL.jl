module OSQPMKL

using OSQP_MKL_jll

export OSQPMKL

"""
    OSQPMKL()

Use the double precision MKL linear algebra backend and linear system solvers.
"""
struct OSQPMKL <: OSQP.OSQPAlgebra end


# The integer type is defined by the C library based on the size of ints on the
# platform being used, so we do the same here.
if Sys.WORD_SIZE == 64   # 64bit system
    const Cc_int = Clonglong
else  # 32bit system
    const Cc_int = Cint
end

###############################################
# Define the library functions                #
###############################################

# Helper functions to wrap the library handles
library(::Type{OSQPMKL}) = OSQP_MKL_jll.osqp_mkl
librarysym(::Type{OSQPMKL}) = :(OSQP_MKL_jll.osqp_mkl)

# Create a definition for all the functions present in the double precision library
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_capabilities
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_version
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_error_message
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_get_dimensions
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_set_default_settings
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_setup
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_solve
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_cleanup
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_warm_start
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_cold_start
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_update_data_vec
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_update_data_mat
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_update_settings
@cdefinition OSQPMKL{Float64} Float64 Cc_int osqp_update_rho

end # module
