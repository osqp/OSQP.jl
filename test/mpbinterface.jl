using MathProgBase
import MathProgBase: LinearQuadraticModel, loadproblem!, setquadobj!, optimize!, status, getobjval, getsolution

@testset "MathProgBase" begin
    solver = OSQPMathProgBaseInterface.OSQPSolver(eps_abs = 1e-7, polish = true)
    MathProgBase.setparameters!(solver, Silent=true)

    # from joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"):
    sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf,solver)
    @test sol.status == :Optimal
    @test isapprox(sol.objval, 130/70, atol=1e-6)
    @test isapprox(norm(sol.sol[1:3] - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)

    @testset "QP1" begin
        m = LinearQuadraticModel(solver)
        loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)

        setquadobj!(m,diagm([10.0,10.0,10.0]))
        rows = [1, 2, 2, 2, 3, 3, 3]
        cols = [1, 1, 1, 2, 2, 3, 3]
        vals = Float64[2, 0.5, 0.5, 2, 1, 1, 1]
        setquadobj!(m,rows,cols,vals)
        optimize!(m)
        stat = status(m)
        @test stat == :Optimal
        @test isapprox(getobjval(m), 130/70, atol=1e-6)
        @test isapprox(norm(getsolution(m) - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)
    end
end

# include(joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"))
# quadprogtest(solver)
# qpdualtest(solver)
