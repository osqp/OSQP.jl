using MathProgBase
import MathProgBase: LinearQuadraticModel, loadproblem!, setquadobj!, optimize!, status, numvar, numconstr, setwarmstart!, getsense,
    getobjval, getsolution, setsense!, getvarLB, setvarLB!, getvarUB, setvarUB!, getconstrLB, setconstrLB!, getconstrUB, setconstrUB!,
    getconstrmatrix, numlinconstr, getsolvetime, getrawsolver, getvartype, setvartype!

@testset "MathProgBase" begin
    solver = OSQPMathProgBaseInterface.OSQPSolver(eps_abs = 1e-7, eps_rel = 1e-16)
    MathProgBase.setparameters!(solver, Silent=true)

    # modified version of joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"):
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
        m2 = copy(m)

        verify_solution = function (m)
            stat = status(m)
            @test stat == :Optimal
            @test isapprox(getobjval(m), 130/70, atol=1e-6)
            @test isapprox(norm(getsolution(m) - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)
            @test getsolvetime(m) > 0
        end

        optimize!(m)
        verify_solution(m)

        setwarmstart!(m2, getsolution(m) .+ 0.1)
        optimize!(m2)
        verify_solution(m2)
    end

    @testset "Unsupported behavior" begin
        m = LinearQuadraticModel(solver)
        @test_throws ErrorException setsense!(m, :Max)
        @test_throws ErrorException loadproblem!(m, spzeros(0, 2), [-Inf, -Inf], [Inf, Inf], [1.; 1.], Float64[], Float64[], :Max) # maximization not supported
        @test_throws ErrorException loadproblem!(m, spzeros(0, 2), [-1., -Inf], [Inf, Inf], [1.; 1.], Float64[], Float64[], :Min) # variable bounds not supported
        @test_throws ErrorException loadproblem!(m, spzeros(0, 2), [-Inf, -Inf], [10., Inf], [1.; 1.], Float64[], Float64[], :Min) # variable bounds not supported

        loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)
        @test_throws ErrorException setvarLB!(m, rand(numvar(m)))
        @test_throws ErrorException setvarUB!(m, rand(numvar(m)))

        @test_throws ErrorException setvartype!(m, [:Cont, :Bin, :Cont])
    end

    @testset "Getters/setters" begin
        m = LinearQuadraticModel(solver)
        loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)

        @test getrawsolver(m) isa OSQP.Model

        @test getsense(m) == :Min

        @test getvarLB(m) == [-Inf, -Inf, -Inf]
        @test getvarUB(m) == [Inf, Inf, Inf]

        @test getconstrLB(m) == [4., 1.]
        setconstrLB!(m, [-4., Inf])
        @test getconstrLB(m) == [-4., Inf]
        @test getconstrUB(m) == [Inf, Inf]
        setconstrUB!(m, [5., 8.])
        @test getconstrUB(m) == [5., 8.]

        @test getconstrmatrix(m) == [1. 2. 3.; 1. 1. 0.]
        @test numlinconstr(m) == 2

        @test getvartype(m) == [:Cont, :Cont, :Cont]
    end
end
