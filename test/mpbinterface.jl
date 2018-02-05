include(joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"))
solver = OSQPMathProgBaseInterface.OSQPSolver(eps_abs = 1e-7, polish = true)
MathProgBase.setparameters!(solver, Silent=true)

quadprogtest(solver)
qpdualtest(solver)
