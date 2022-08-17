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
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_capabilities
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_version
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_error_message
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_get_dimensions
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_set_default_settings
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_setup
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_solve
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_cleanup
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_warm_start
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_cold_start
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_update_data_vec
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_update_data_mat
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_update_settings
@cdefinition OSQPCUDAAlgebra{Float64} osqp_cuda_double Float64 Cint osqp_update_rho

# Create a definition for all the functions present in the single precision library
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_capabilities
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_version
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_error_message
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_get_dimensions
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_set_default_settings
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_setup
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_solve
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_cleanup
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_warm_start
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_cold_start
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_update_data_vec
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_update_data_mat
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_update_settings
@cdefinition OSQPCUDAAlgebra{Float32} osqp_cuda_single Float32 Cint osqp_update_rho


end # module
