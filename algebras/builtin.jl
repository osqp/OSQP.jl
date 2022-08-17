using OSQP_jll

export OSQPBuiltinAlgebra

# The integer type is defined by the C library based on the size of ints on the
# platform being used, so we do the same here.
if Sys.WORD_SIZE == 64   # 64bit system
    const Cc_int = Clonglong
else  # 32bit system
    const Cc_int = Cint
end

"""
    OSQPBuiltinAlgebra{FT<:Union{Float32, Float64}} <: OSQPAlgebra{FT, Cc_int}

Use the builtin OSQP linear algebra functions and linear system solvers.

This algebra backend supports either single or double precision, which can
be chosen by specifying either `Float32` or `Float64`, respectively, as the
type parameter.
"""
struct OSQPBuiltinAlgebra{FT<:Union{Float32, Float64}} <: OSQPAlgebra{FT, Cc_int}
end

"""
    OSQPBuiltinAlgebra()

Use the builtin double precision linear algebra functions and linear system solvers.
"""
function OSQPBuiltinAlgebra()
    return OSQPBuiltinAlgebra{Float64}()
end

###############################################
# Define the double precision library         #
###############################################
# Create a definition for all the functions present in the double precision library
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_capabilities
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_version
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_error_message
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_get_dimensions
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_set_default_settings
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_setup
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_solve
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_cleanup
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_warm_start
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_cold_start
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_update_data_vec
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_update_data_mat
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_update_settings
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double Float64 Cc_int osqp_update_rho



###############################################
# Define the single precision library         #
###############################################
# Create a definition for all the functions present in the double precision library
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_capabilities
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_version
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_error_message
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_get_dimensions
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_set_default_settings
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_setup
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_solve
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_cleanup
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_warm_start
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_cold_start
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_update_data_vec
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_update_data_mat
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_update_settings
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single Float32 Cc_int osqp_update_rho
