function setup_feasibility()
    options = Dict(:verbose => false,
                   :eps_abs => 1e-06,
                   :eps_rel => 1e-06,
                   :max_iter => 5000)
    return options
end

tol = 1e-3


@testset "feasibility" begin

    @testset "feasibility_problem" begin
        n = 30
        m = 30
        P = spzeros(n, n)
        q = zeros(n)
        A = sprandn(m, n, 0.8)
        u = randn(m)
        l = u
        options = setup_feasibility()

        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)

        @test isapprox(norm(A * results.x - u), 0., atol=tol)


    end

end
