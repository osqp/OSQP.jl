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
Tdoubledict = Dict(
    :Tfloat => :Float64,
    :Tint => :Cc_int
)
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_capabilities          Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_version               Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_error_message         Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_get_dimensions        Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_set_default_settings  Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_setup                 Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_solve                 Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_cleanup               Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_warm_start            Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_cold_start            Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_update_data_vec       Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_update_data_mat       Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_update_settings       Tdoubledict
@cdefinition OSQPMKLAlgebra osqp_mkl osqp_update_rho            Tdoubledict

end # module
