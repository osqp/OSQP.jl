using OSQP
using Compat
using Compat.Test, Compat.SparseArrays, Compat.LinearAlgebra, Compat.Random

if isdefined(Random, :seed!)
    using Random: seed!
else
    const seed! = srand
end

tests = [
    "basic.jl",
    "dual_infeasibility.jl",
    "feasibility.jl",
    "polishing.jl",
    "primal_infeasibility.jl",
    "unconstrained.jl",
    "warm_start.jl",
    "update_matrices.jl",
    "mpbinterface.jl",
    "MathOptInterfaceOSQP.jl"
    ]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end