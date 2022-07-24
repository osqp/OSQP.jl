using OSQP_jll

export OSQPBuiltin

"""
    OSQPBuiltin{FT} <: OSQPAlgebra where {FT <: Union{Float32, Float64}}

Use the builtin OSQP linear algebra functions and linear system solvers.

This algebra backend supports either single or double precision, which can
be chosen by specifying either `Float32` or `Float64`, respectively, as the
type parameter.
"""
struct OSQPBuiltin{FT} <: OSQPAlgebra where {FT <: Union{Float32, Float64}}
end

"""
    OSQPBuiltin()

Use the builtin double precision linear algebra functions and linear system solvers.
"""
function OSQPBuiltin()
    return OSQPBuiltin{Float64}()
end

# The integer type is defined by the C library based on the size of ints on the
# platform being used, so we do the same here.
if Sys.WORD_SIZE == 64   # 64bit system
    const Cc_int = Clonglong
else  # 32bit system
    const Cc_int = Cint
end

###############################################
# Define the double precision library         #
###############################################

# Helper functions to wrap the library handles
library(::Type{OSQPBuiltin{Float64}}) = OSQP_jll.osqp_builtin_double
librarysym(::Type{OSQPBuiltin{Float64}}) = :(OSQP_jll.osqp_builtin_double)

# Create a definition for all the functions present in the double precision library
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_capabilities
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_version
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_error_message
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_get_dimensions
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_set_default_settings
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_setup
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_solve
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_cleanup
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_warm_start
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_cold_start
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_update_data_vec
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_update_data_mat
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_update_settings
@cdefinition OSQPBuiltin{Float64} Float64 Cc_int osqp_update_rho



###############################################
# Define the single precision library         #
###############################################

# Helper functions to wrap the library handles
library(::Type{OSQPBuiltin{Float32}}) = OSQP_jll.osqp_builtin_single
librarysym(::Type{OSQPBuiltin{Float32}}) = :(OSQP_jll.osqp_builtin_single)

# Create a definition for all the functions present in the double precision library
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_capabilities
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_version
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_error_message
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_get_dimensions
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_set_default_settings
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_setup
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_solve
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_cleanup
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_warm_start
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_cold_start
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_update_data_vec
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_update_data_mat
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_update_settings
@cdefinition OSQPBuiltin{Float32} Float32 Cc_int osqp_update_rho
