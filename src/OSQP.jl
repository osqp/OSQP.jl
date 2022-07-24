module OSQP

export OSQPMathProgBaseInterface
using SparseArrays
using LinearAlgebra


include("constants.jl")
include("types.jl")

include("algebra_utils.jl")
include("prototypes.jl")

# This is the only algebra that is included by default in the OSQP package
# The other algebras are provided as other packages
include("../algebras/builtin.jl")

include("interface.jl")
include("MOI_wrapper.jl")
const Optimizer = MathOptInterfaceOSQP.Optimizer

end # module
