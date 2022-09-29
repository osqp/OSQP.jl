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
Tdoubledict = Dict(
    :Tfloat => :Float64,
    :Tint => :Cc_int
)
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_capabilities          Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_version               Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_error_message         Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_get_dimensions        Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_set_default_settings  Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_setup                 Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_solve                 Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_cleanup               Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_warm_start            Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_cold_start            Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_update_data_vec       Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_update_data_mat       Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_update_settings       Tdoubledict
@cdefinition OSQPBuiltinAlgebra{Float64} osqp_builtin_double osqp_update_rho            Tdoubledict



###############################################
# Define the single precision library         #
###############################################
Tsingledict = Dict(
    :Tfloat => :Float32,
    :Tint => :Cc_int
)
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_capabilities          Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_version               Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_error_message         Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_get_dimensions        Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_set_default_settings  Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_setup                 Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_solve                 Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_cleanup               Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_warm_start            Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_cold_start            Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_update_data_vec       Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_update_data_mat       Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_update_settings       Tsingledict
@cdefinition OSQPBuiltinAlgebra{Float32} osqp_builtin_single osqp_update_rho            Tsingledict
