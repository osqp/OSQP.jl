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
        using Random
        rng = Random.MersenneTwister(666)


        n = 30
        m = 50
        Pt = sprandn(rng, n, n, 0.5)
        P = Pt * Pt'
        q = randn(rng, n)
        A = sprandn(rng, m, n, 0.5)
        l = -3 .+ randn(rng, m)
        u = 3 .+ randn(rng, m)

        # Solve problem
        options = setup_polishing()
        model = OSQP.Model()
        OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
        results = OSQP.solve!(model)


        # # Solve with Mosek for reference solution
        # using JuMP, MosekTools
        #
        # jump_model = JuMP.Model(with_optimizer(Mosek.Optimizer));
        # @variable(jump_model, x[1:n]);
        # @objective(jump_model, Min, 0.5 * x' * P * x  + q' * x);
        # con = @constraint(jump_model, [i = 1:m], l[i] <= Vector(A[i, :])' * x <= u[i]);
        #
        # status = JuMP.optimize!(jump_model);
        # obj_test = JuMP.objective_value(jump_model);
        # x_test = JuMP.value.(x);
        # y_test = -JuMP.dual.(con)
        # println(x_test)
        # println(y_test)
        # println(obj_test)

        x_test = [0.41137084480890623, 0.28318004895183285, 0.4777528132318454, -0.23164842615545447, -0.3629757238797728, 0.14431287989072183, 0.6438564504118605, -0.5761094276277431, -0.05982787893150773, 0.15030576101377005, 0.028787631121672347, -0.27216082894024013, 0.17141240340176417, -0.273743108255697, 0.15630235714216442, -0.06763333754014154, -1.020331708945701, -0.012617123323538763, -0.6683845746964124, 0.5531428372752116, -0.1320489823459991, -0.36749007989638366, 0.0595756204214971, -0.2428072482922348, 0.6357746681280918, 0.9216345137918707, -0.14249830279317166, 0.3440205926304816, -0.2570648319065089, -0.6783094533788961];


        y_test = [-1.0282863481250939e-13, -3.37258508945909e-14, 3.4296811845328215e-14, 1.777147102264673e-14, -2.4893358571930445e-14, 1.322599809069287e-14, -1.6595755569514817e-13, -8.317794666110387e-14, -9.2640088352221e-15, 1.253500392114895e-13, 1.0854356507664321e-13, 3.1706833729517506e-14, -6.241276348174211e-13, -8.535981813127953e-14, 3.5093009175770576e-13, 7.29118035417498e-14, 3.618176044780794e-14, -5.491504492855572e-15, -5.667489500506908e-14, 0.0050306281848406156, 5.067317066489025e-14, -1.0803519472826196e-14, 5.237647375580867e-14, 9.474300815751858e-14, 2.9853746376673826e-14, -1.1895032769230835e-14, -1.3441509011114794e-13, -0.08678665919072658, -0.04228028938433161, -6.59489689104704e-14, -5.710323913717883e-14, 5.390342831441443e-14, -6.35787374994456e-14, 1.252381079606396e-13, 0.22585055797097334, -0.9388516056549417, 1.0170683219018741e-13, 0.2625913401764051, -3.454237726808208e-13, -0.08487116932540341, 5.24402777336478e-14, 4.099664389450276e-14, 3.614515769510322e-14, 2.1517424811393735e-14, 1.5159680259859412e-14, -4.016603096658772e-14, 5.0318209747657e-13, -2.981549776801216e-14, -1.1248637104537303e-14, -3.8425095057919334e-13]

        obj_test = -6.0748863556792045

        @test isapprox(results.x, x_test, atol=tol)
        @test isapprox(results.y, y_test, atol=tol)
        @test isapprox(results.info.obj_val, obj_test, atol=tol)
        @test results.info.status_polish == 1

    end
end
