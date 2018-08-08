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
        P = sparse(Diagonal(randn(n)) + 1. * sparse(I, n, n))
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
        seed!(1)

        n = 30
        m = 50
        Pt = sprandn(n, n, 0.5)
        P = Pt * Pt'
        q = randn(n)
        A = sprandn(m, n, 0.5)
        l = -3 .+ randn(m)
        u = 3 .+ randn(m)

        # Solve problem
        options = setup_polishing()
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)


        # # Solve with Gurobi
        # using Gurobi
        # env = Gurobi.Env()
        # setparam!(env, "OutputFlag", 1)
        # model = Gurobi.Model(env, "model")
        # add_cvars!(model, q); update_model!(model)
        # add_constrs!(model, A, repmat(['<'], m), u); update_model!(model)
        # add_constrs!(model, A, repmat(['>'], m), l); update_model!(model)
        # add_qpterms!(model, P); update_model!(model)
        # optimize(model)
        # x_test = get_solution(model)
        # obj_test= get_objval(model)
        # y_test = -Gurobi.get_dblattrarray(model, "Pi", 1, 2 * m)
        # y_test = y_test[m+1:end] + y_test[1:m]
                #
        # println(x_test)
        # println(y_test)
        # println(obj_test)

        x_test = [0.369834, -0.0209277, 0.068939, -0.604151, 0.60773, -0.715965, -0.128837, -0.593642, 0.928932, 0.678794, 0.194334, 0.228124, 0.188779, 0.607246, 0.479752, 0.0878951, 0.275489, 0.761593, 0.0600479, -0.950663, -0.191354, 1.12286, -0.593261, -0.0932879, 0.0286826, -0.0257784, 0.358456, -0.0693514, -0.00358498, -0.0270959]
        y_test = [4.88779e-14, 0.441478, 5.89704e-14, 1.76439e-13, 8.10529e-15, 4.9186e-14, 5.00637e-14, 2.00889e-14, 5.39402e-13, -1.9796e-14, 1.52041e-13, 0.149979, 8.97844e-12, 5.93941e-14, -4.11719e-14, -2.38358e-13, 0.226114, -1.8933e-14, 1.54839e-13, 2.99422e-14, 7.16687e-15, -5.32536e-14, 0.116129, -2.53245e-14, 1.21586e-14, -1.66761e-13, 3.04905e-14, 7.35909e-14, 3.84878e-12, -3.09855e-15, 4.30889e-9, 4.09921e-14, -2.4863e-14, -2.46841e-14, 4.37156e-15, 0.219029, -5.14155e-14, 8.31705e-14, 4.51252e-15, -7.90941e-14, 0.206435, -0.196264, 1.06762e-14, 0.338783, 7.58781e-10, -0.105739, 3.30014e-13, -4.05467e-14, 3.3212e-14, -0.161755]

        obj_test = -6.021607803961248

        @test isapprox(results.x, x_test, atol=tol)
        @test isapprox(results.y, y_test, atol=tol)
        @test isapprox(results.info.obj_val, obj_test, atol=tol)
        @test results.info.status_polish == 1

    end
end
