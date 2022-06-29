function setup_primal_infeasibility()
    options = Dict(
        :verbose => false,
        :eps_abs => 1e-05,
        :eps_rel => 1e-05,
        :eps_dual_inf => 1e-18,
        :scaling => true,
    )
    return options
end

tol = 1e-5

@testset "primal_infeasibility" begin
    @testset "primal_infeasible_problem" begin
        seed!(1)

        n = 50
        m = 500
        P = sprandn(n, n, 0.6)
        P = P' * P
        q = randn(n)
        A = sprandn(m, n, 0.6)
        u = 3 .+ randn(m)
        l = -3 .+ randn(m)

        # Make problem infeasible
        A[Int(n / 2), :] = A[Int(n / 2)+1, :]
        l[Int(n / 2)] = u[Int(n / 2)+1] + 10 * rand()
        u[Int(n / 2)] = l[Int(n / 2)] + 0.5

        options = setup_primal_infeasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P = P, q = q, A = A, l = l, u = u, options...)
        results = OSQP.solve!(model)

        @test results.info.status == :Primal_infeasible
    end

    @testset "primal_dual_infeasible_problem" begin
        seed!(1)

        n = 2
        m = 4
        P = spzeros(n, n)
        q = [-1.0, -1]
        A = sparse([1.0 -1; -1.0 1.0; 1.0 0.0; 0.0 1.0])
        l = [1.0, 1.0, 0.0, 0.0]
        u = Inf * ones(m)

        options = setup_primal_infeasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P = P, q = q, A = A, l = l, u = u, options...)
        results = OSQP.solve!(model)

        @test results.info.status == :Primal_infeasible
    end
end
