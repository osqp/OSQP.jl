using OSQP
using Base.Test


# tests = ["simple.jl"]
tests = ["basic_tests.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end


# TODO: Add QP tests MathProgBase
