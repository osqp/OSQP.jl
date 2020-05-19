function setup_polishing()
    options = Dict(:verbose => false,
                   :polish => true,
                   :eps_abs => 1e-03,
                   :eps_rel => 1e-03,
                   :verbose => false,
                   :max_iter => 5000)
    return options
end

tol = 1e-3


@testset "polishing" begin

    @testset "polishing_problem" begin
        P = sparse(Diagonal([11.; 0.]))
        q = [3.; 4.]
        A = sparse([-1. 0.; 0. -1.; -1. -3; 2. 5.; 3. 4.])
        u = [0.; 0.; -15.; 100.; 80]
        l = -Inf * ones(length(u))
        (m, n) = size(A)

        # Solve problem
        options = setup_polishing()
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)


        x_test = [9.90341e-11; 5.0]
        y_test = [1.66667; 0.0; 1.33333; 1.20431e-14; 1.49741e-14]
        obj_test = 20.

        @test isapprox(results.x, x_test, atol=tol)
        @test isapprox(results.y, y_test, atol=tol)
        @test isapprox(results.info.obj_val, obj_test, atol=tol)
        @test results.info.status_polish == 1


    end


    @testset "polishing_unconstrained" begin

        seed!(1)

        n = 10
        m = n
        P = sparse(Diagonal(rand(n)) + .2 * sparse(I, n, n))
        q = randn(n)
        A = sparse(I, n, n)
        l = -100 * ones(m)
        u = 100 * ones(m)

        # Solve problem
        options = setup_polishing()
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)


        # Explicit solution
        invP = inv(Array(P))
        x_test = - invP * q
        y_test = zeros(m)
        obj_test = - .5 * q' * invP * q

        @test isapprox(results.x, x_test, atol=tol)
        @test isapprox(results.y, zeros(m), atol=tol)
        @test isapprox(results.info.obj_val, obj_test, atol=tol)
        @test results.info.status_polish == 1
    end

    @testset "polish_random" begin

        # load randomly generated problem with known accurate solution from Mosek
        problem_data = FileIO.load("./problem_data/random_polish_qp.jld2");
        P = problem_data["P"]; q = problem_data["q"]; A = problem_data["A"]; u = problem_data["u"]; l = problem_data["l"];
        options = setup_polishing()
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)


        @test isapprox(results.x, problem_data["x_test"], atol=tol)
        @test isapprox(results.y, problem_data["y_test"], atol=tol)
        @test isapprox(results.info.obj_val, problem_data["obj_test"], atol=tol)
        @test results.info.status_polish == 1

    end
end
