module OSQPCUDA

using OSQP_CUDA_jll

using OSQP
using OSQP: OSQPAlgebra
using OSQP: @cdefinition

export OSQPCUDAAlgebra

"""
    OSQPCUDAAlgebra{FT} <: OSQPAlgebra where {FT <: Union{Float32, Float64}}

Use the CUDA OSQP linear algebra functions and linear system solvers.

This algebra backend supports either single or double precision, which can
be chosen by specifying either `Float32` or `Float64`, respectively, as the
type parameter.
"""
struct OSQPCUDAAlgebra{FT <: Union{Float32, Float64}} <: OSQPAlgebra{FT,Cint}
end

"""
    OSQPCUDAAlgebra()

Use the CUDA double precision linear algebra functions and linear system solvers.
"""
function OSQPCUDAAlgebra()
    return OSQPCUDAAlgebra{Float64}()
end

###############################################
# Define the library functions                #
###############################################
# Create a definition for all the functions present in the double precision library
Tdoubledict = Dict(
    :Tfloat => :Float64,
    :Tint => :Cint
)
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_capabilities          Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_version               Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_error_message         Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_get_dimensions        Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_set_default_settings  Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_setup                 Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_solve                 Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_cleanup               Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_warm_start            Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_cold_start            Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_update_data_vec       Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_update_data_mat       Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_update_settings       Tdoubledict
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double osqp_update_rho            Tdoubledict


# Create a definition for all the functions present in the single precision library
Tsingledict = Dict(
    :Tfloat => :Float32,
    :Tint => :Cint
)
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_capabilities          Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_version               Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_error_message         Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_get_dimensions        Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_set_default_settings  Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_setup                 Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_solve                 Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_cleanup               Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_warm_start            Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_cold_start            Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_update_data_vec       Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_update_data_mat       Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_update_settings       Tsingledict
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single osqp_update_rho            Tsingledict


end # module
