using OSQP
using Base.Test


tests = ["basic.jl", 
	 "dual_infeasibility.jl",
	 "feasibility.jl",
	 "polishing.jl",
	 "primal_infeasibility.jl",
	 "unconstrained.jl",
	 "warm_start.jl",
	 "update_matrices.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end


# TODO: Add QP tests MathProgBase
