module OSQPMKL

using OSQP_MKL_jll

using OSQP
using OSQP: OSQPAlgebra
using OSQP: @cdefinition

export OSQPMKLAlgebra

"""
    OSQPMKLAlgebra()

Use the double precision MKL linear algebra backend and linear system solvers.
"""
struct OSQPMKLAlgebra <: OSQPAlgebra end


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
# Create a definition for all the functions present in the double precision library
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_capabilities
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_version
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_error_message
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_get_dimensions
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_set_default_settings
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_setup
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_solve
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_cleanup
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_warm_start
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_cold_start
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_update_data_vec
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_update_data_mat
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_update_settings
@cdefinition OSQPMKLAlgebra osqp_mkl Float64 Cc_int osqp_update_rho

end # module
