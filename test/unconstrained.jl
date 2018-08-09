function setup_unconstrained()
    options = Dict(:verbose => false,
                   :eps_abs => 1e-08,
                   :eps_rel => 1e-08,
                   :eps_dual_inf => 1e-18)
    return options
end

tol = 1e-5


@testset "unconstrained" begin

    @testset "unconstrained_problem" begin
        seed!(1)

        n = 30
        m = 0
        P = sparse(Diagonal(rand(n)) + 0.2 * sparse(I, n, n))
        q = randn(n)
        A = spzeros(m, n)
        u = Float64[]
        l = Float64[]


        # Explicit solution
        invP = inv(Array(P))
        x_test = - invP * q
        y_test = zeros(m)
        obj_test = - .5 * q' * invP * q

        options = setup_unconstrained()

        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)

        @test isapprox(results.x, x_test, atol=tol)
        @test isapprox(results.y, y_test, atol=tol)
        @test isapprox(results.info.obj_val, obj_test, atol=tol)
        @test results.info.status == :Solved

    end

end
