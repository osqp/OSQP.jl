using MathProgBase
using MathProgBase.SolverInterface


# Include linprog.jl tests
if isdefined(Base, :pathof)
    include(joinpath(dirname(pathof(MathProgBase)), "..", "test", "linprog.jl"))
    include(joinpath(dirname(pathof(MathProgBase)), "..", "test", "linproginterface.jl"))
else
    include(joinpath(Pkg.dir("MathProgBase"), "test", "linprog.jl"))
    include(joinpath(Pkg.dir("MathProgBase"), "test", "linproginterface.jl"))
end

@testset "MathProgBase" begin
    solver = OSQPMathProgBaseInterface.OSQPSolver(eps_abs = 1e-7, eps_rel = 1e-16)
    MathProgBase.setparameters!(solver, Silent=true)

    @testset "linprog" begin
        linprogtest(solver)
    end

    @testset "linproginterface" begin
        linprogsolvertest(solver, 1e-03)
    end

    @testset "quadprog" begin
        # modified from joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"):
        sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf,solver)
        @test sol.status == :Optimal
        @test isapprox(sol.objval, 130/70, atol=1e-6)
        @test isapprox(norm(sol.sol[1:3] - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)
    end

    @testset "QP1" begin
        # modified from joinpath(Pkg.dir("MathProgBase"), "test", "quadprog.jl"):
        m = LinearQuadraticModel(solver)
        loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)

        setquadobj!(m,diagm(0 => [10.0,10.0,10.0]))
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

    @testset "basic_QP" begin
        # same QP as test/basic.jl: "basic_QP" but using MathProgBase
        P = sparse([11. 0.; 0. 0.])
        q = [3.; 4]
        A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4])
        u = [0.; 0.; -15; 100; 80]
        l = fill(-Inf, length(u))

        verify_solution = function (m)
            tol = 1e-5

            @test isapprox(norm(getsolution(m) - [0.; 5.]), 0., atol=tol)
            @test isapprox(norm(getconstrduals(m) - [-1.666666666666; 0.; -1.3333333; 0.; 0.]), 0., atol=tol)
            @test isapprox(getobjval(m), 20., atol=tol)
        end

        m1 = LinearQuadraticModel(solver)
        loadproblem!(m1, A, [-Inf, -Inf], [Inf, Inf], q, l, u, :Min)
        m2 = copy(m1)
        setquadobj!(m1, P)
        optimize!(m1)
        verify_solution(m1)

        setquadobj!(m2, findnz(triu(P))...) # triu to avoid duplicate elements
        optimize!(m2)
        verify_solution(m2)
    end

    @testset "Unsupported behavior" begin
        m = LinearQuadraticModel(solver)
        loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)
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
