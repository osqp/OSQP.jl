
@testset "non_convex" begin

    @testset "non_convex_small_sigma" begin

        # Setup problem
        P = sparse([2. 5.; 5. 1.])
        q = [3.; 4]
        A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4])
        u = [0.; 0.; -15; 100; 80]
        l = -Inf * ones(length(u))
        options = Dict(:verbose => false, :sigma => 1e-06)
        model = OSQP.Model()
        try
            # Setup should fail due to (P + sigma I) having a negative eigenvalue
            global test_setup = 1
            OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        catch
            global test_setup = 0
        end

        @test test_setup == 0

    end


    @testset "non_convex_big_sigma" begin

        # Setup workspace with new sigma
        P = sparse([2. 5.; 5. 1.])
        q = [3.; 4]
        A = sparse([-1 0; 0 -1; -1 -3; 2 5; 3 4])
        u = [0.; 0.; -15; 100; 80]
        l = -Inf * ones(length(u))
        options = Dict(:verbose => false, :sigma => 5.)
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)

        # Solve problem
        results = OSQP.solve!(model)

        @test isnan(results.info.obj_val)
        @test results.info.status == :Non_convex

    end

end
