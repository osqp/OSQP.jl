using OSQP

s = OSQPSolver()



using JuMP
# MathProgBase.setparameters!(s; Silent = true)
m = Model(solver=s)
@variable(m, x)
@variable(m, y)
@constraint(m, 1 <= x <= 4)
@constraint(m, y >= 3)
@objective(m, Min, x^2 + y^2)
status = solve(m)
@show getvalue(x)
@show getvalue(y)
@show getobjectivevalue(m)



# MathProg tests as in Gurobi.jl
# include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
# linprogtest(s)

# include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
# linprogsolvertest(s)

# include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
# quadprogtest(s)
# qpdualtest(s)
