function setup_dual_infeasibility()
        options = Dict(:verbose => false,
                       :eps_abs => 1e-05,
                       :eps_rel => 1e-05,
               :eps_prim_inf => 1e-15,
               :check_termination => 1)
    return options
end

tol = 1e-5


@testset "dual_infeasibility" begin

    @testset "dual_infeasible_lp" begin
        P = spzeros(2, 2)
        q = [2.; -1.]
        A = sparse(I, 2, 2)
        u = Inf * ones(2)
        l = [0.; 0.]
        options = setup_dual_infeasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)

        @test results.info.status == :Dual_infeasible

    end

    @testset "dual_infeasible_qp" begin
        P = sparse(Diagonal([4.; 0.]))
        q = [0.; 2]
        A = sparse([1. 1.; -1. 1])
        u = [2.; 3]
        l = -Inf * ones(2)
        options = setup_dual_infeasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)

        @test results.info.status == :Dual_infeasible

    end

    @testset "primal_dual_infeasible" begin
        P = spzeros(2, 2)
        q = [-1.; -1.]
        A = sparse([1. -1.; -1. 1; 1. 0; 0. 1])
        u = Inf * ones(4)
        l = [1., 1., 0., 0.]
        options = setup_dual_infeasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)

        # Warm start to avoid infeasibility detection at first step
        OSQP.warm_start!(model, x=[50.; 30.], y=[-2.; -2; -2; -2])

        results = OSQP.solve!(model)

        @test results.info.status == :Dual_infeasible

    end
end
