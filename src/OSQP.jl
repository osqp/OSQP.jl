module OSQP

export OSQPMathProgBaseInterface
using SparseArrays
using LinearAlgebra

using OSQP_jll

include("constants.jl")
include("types.jl")
include("interface.jl")
include("MOI_wrapper.jl")
const Optimizer = MathOptInterfaceOSQP.Optimizer

end # module
