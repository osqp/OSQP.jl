using OSQP
using Base.Test


tests = ["basic.jl", 
	 "dual_infeasibility.jl",
	 "feasibility.jl",
	 "polishing.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end


# TODO: Add QP tests MathProgBase
