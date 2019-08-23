using OSQP
using OSQP.MathOptInterfaceOSQP
using LinearAlgebra
using Random
using SparseArrays
using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

const Affine = MOI.ScalarAffineFunction{Float64}

@testset "ProblemModificationCache" begin
    rng = MersenneTwister(1234)
    n = 15
    m = 10
    q = randn(rng, n)
    P = (X = sprandn(rng, n, n, 0.1); X' * X)
    P += eps() * I # needed for a test later on
    l = -rand(rng, m)
    u = rand(rng, m)
    A = sprandn(rng, m, n, 0.6)
    modcache = MathOptInterfaceOSQP.ProblemModificationCache(P, q, A, l, u)

    model = OSQP.Model()
    OSQP.setup!(model; P = P, q = q, A = A, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    baseresults = OSQP.solve!(model)
    @test baseresults.info.status == :Solved

    # Modify q, ensure that updating results in the same solution as calling setup! with the modified q
    @test !modcache.q.dirty
    modcache.q[3] = 5.
    @test modcache.q.dirty
    @test modcache.q.data[3] == 5.
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    @test !modcache.q.dirty
    qmod_update_results = OSQP.solve!(model)
    @test !isapprox(baseresults.x, qmod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    model2 = OSQP.Model()
    OSQP.setup!(model2; P = P, q = modcache.q.data, A = A, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    qmod_setup_results = OSQP.solve!(model2)
    @test qmod_update_results.x ≈ qmod_setup_results.x atol = 1e-7

    # Modify A, ensure that updating results in the same solution as calling setup! with the modified A and q
    (rows, cols, _) = findnz(A)
    Amodindex = rand(rng, 1 : nnz(A))
    row = rows[Amodindex]
    col = cols[Amodindex]
    val = randn(rng)
    modcache.A[row, col] = val
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    @test isempty(modcache.A.modifications)
    Amod_update_results = OSQP.solve!(model)
    @test !isapprox(baseresults.x, Amod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    @test !isapprox(qmod_update_results.x, Amod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    Amod = copy(A)
    Amod[row, col] = val
    model3 = OSQP.Model()
    OSQP.setup!(model3; P = P, q = modcache.q.data, A = Amod, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    Amod_setup_results = OSQP.solve!(model3)
    @test Amod_update_results.x ≈ Amod_setup_results.x atol = 1e-7

    # MatrixModificationCache: colon indexing
    modcache.P[:] = 0.
    for i = 1 : n
        modcache.P[i, i] = 1.
    end
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    Pmod_update_results = OSQP.solve!(model)
    model4 = OSQP.Model()
    Pmod = sparse(1.0I, n, n)
    OSQP.setup!(model4; P = Pmod, q = modcache.q.data, A = Amod, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    Pmod_setup_results = OSQP.solve!(model4)
    @test Pmod_update_results.x ≈ Pmod_setup_results.x atol = 1e-8

    # Modifying the sparsity pattern is not allowed
    nzinds = map(CartesianIndex, zip(rows, cols))
    zinds = zinds = setdiff(vec(CartesianIndices(A)), nzinds)
    for zind in zinds
        @test_throws ArgumentError modcache.A[zind[1], zind[2]] = randn(rng)
    end
    @test_throws ArgumentError modcache.A[:] = 1
end

# FIXME: type piracy. Generalize and move to MOIU.
function MOI.get(optimizer::MOIU.CachingOptimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction, <:Any})
    f = MOI.get(optimizer, MOI.ConstraintFunction(), ci)
    ret = f.constant
    for term in f.terms
        ret += term.coefficient * MOI.get(optimizer, MOI.VariablePrimal(), term.variable_index)
    end
    ret
end

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

function defaultoptimizer()
    optimizer = OSQP.Optimizer()
    MOI.set(optimizer, MOI.Silent(), true)
    MOI.set(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOI.set(optimizer, OSQPSettings.MaxIter(), 10000)
    MOI.set(optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    return optimizer
end
function bridged_optimizer()
    optimizer = defaultoptimizer()
    cached = MOIU.CachingOptimizer(MOIU.UniversalFallback(OSQPModel{Float64}()), optimizer)
    return MOI.Bridges.full_bridge_optimizer(cached, Float64)
end

@testset "CachingOptimizer: unit" begin
    excludes = [# TODO
                "time_limit_sec",
                # Quadratic constraints are not supported
                "solve_qcp_edge_cases",
                # No method get(::Optimizer, ::MathOptInterface.ConstraintPrimal, ::MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.Nonpositives})
                "solve_duplicate_terms_vector_affine",
                # FIXME KeyError: key CartesianIndex(1, 2) not found
                "solve_qp_edge_cases",
                # ConstraintPrimal not supported
                "solve_affine_deletion_edge_cases",
                # Integer and ZeroOne sets are not supported
                "solve_integer_edge_cases", "solve_objbound_edge_cases",
                "solve_zero_one_with_bounds_1",
                "solve_zero_one_with_bounds_2",
                "solve_zero_one_with_bounds_3"]

    MOIT.unittest(bridged_optimizer(), config, excludes)
end

@testset "CachingOptimizer: linear problems" begin
    excludes = ["partial_start"]  # See comment https://github.com/JuliaOpt/MathOptInterface.jl/blob/ecf691545e67552ff437ed26ec4ddfff03c50327/src/Test/contlinear.jl#L1715
    append!(excludes,
    if Int == Int32
        ["linear7"] # https://github.com/JuliaOpt/MathOptInterface.jl/issues/377#issuecomment-394912761
    else
        []
    end)
    MOIT.contlineartest(bridged_optimizer(), config, excludes)
end

@testset "CachingOptimizer: quadratic problems" begin
    excludes = String[]
    optimizer = defaultoptimizer()
    MOIT.qptest(bridged_optimizer(), config, excludes)
end

function test_optimizer_modification(modfun::Base.Callable, model::MOI.ModelLike, optimizer::T, idxmap::MOIU.IndexMap,
        cleanoptimizer::T, config::MOIT.TestConfig) where T<:MOI.AbstractOptimizer
    # apply modfun to both the model and the optimizer
    modfun(model)
    modfun(optimizer)

    # copy model into clean optimizer
    cleanidxmap = MOI.copy_to(cleanoptimizer, model)
    @test cleanidxmap.varmap == idxmap.varmap
    @test cleanidxmap.conmap == idxmap.conmap

    # call optimize! on both optimizers
    MOI.optimize!(optimizer)
    MOI.optimize!(cleanoptimizer)

    # compare results
    atol = config.atol
    rtol = config.rtol
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.get(cleanoptimizer, MOI.TerminationStatus())
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.get(cleanoptimizer, MOI.PrimalStatus())
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ MOI.get(cleanoptimizer, MOI.ObjectiveValue()) atol=atol rtol=rtol

    if MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        modelvars = MOI.get(model, MOI.ListOfVariableIndices())
        for v_model in modelvars
            v_optimizer = idxmap[v_model]
            @test MOI.get(optimizer, MOI.VariablePrimal(), v_optimizer) ≈ MOI.get(cleanoptimizer, MOI.VariablePrimal(), v_optimizer) atol=atol rtol=rtol
        end
    end

    if config.duals
        @test MOI.get(optimizer, MOI.DualStatus()) == MOI.get(cleanoptimizer, MOI.DualStatus())
        if MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
            for (F, S) in MOI.get(model, MOI.ListOfConstraints())
                cis_model = MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
                for ci_model in cis_model
                    ci_optimizer = idxmap[ci_model]
                    @test MOI.get(optimizer, MOI.ConstraintDual(), ci_optimizer) ≈ MOI.get(cleanoptimizer, MOI.ConstraintDual(), ci_optimizer) atol=atol rtol=rtol
                end
            end
        end
    end
end

function zero_warm_start!(optimizer::MOI.ModelLike, vars, cons)
    for vi in vars
        MOI.set(optimizer, MOI.VariablePrimalStart(), vi, 0.0)
    end
    for ci in cons
        MOI.set(optimizer, MOI.ConstraintDualStart(), ci, -0.0)
    end
end

term(c, x::MOI.VariableIndex) = MOI.ScalarAffineTerm(c, x)
term(c, x::MOI.VariableIndex, y::MOI.VariableIndex) = MOI.ScalarQuadraticTerm(c, x, y)

@testset "No CachingOptimizer: problem modification after copy_to" begin
    # Initial setup: modified version of MOIT.linear1test
    # min -x
    # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
    #       x, y >= 0   (x, y ∈ Nonnegatives)
    model = OSQPModel{Float64}()
    MOI.empty!(model)
    v = MOI.add_variables(model, 2)
    x, y = v
    cf = MOI.ScalarAffineFunction([term.([0.0, 0.0], v); term.([1.0, 1.0], v); term.([0.0, 0.0], v)], 0.0)
    c = MOI.add_constraint(model, cf, MOI.Interval(-Inf, 1.0))
    vc1 = MOI.add_constraint(model, MOI.SingleVariable(v[1]), MOI.Interval(0.0, Inf))
    vc2 = MOI.add_constraint(model, v[2], MOI.Interval(0.0, Inf))
    objf = MOI.ScalarAffineFunction([term.([0.0, 0.0], v); term.([-1.0, 0.0], v); term.([0.0, 0.0], v)], 0.0)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objf)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    optimizer = defaultoptimizer()
    idxmap = MOI.copy_to(optimizer, model)
    @test MOI.get(optimizer, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 2
    @test MOI.get(optimizer, MOI.ListOfVariableIndices()) == [MOI.VariableIndex(1), MOI.VariableIndex(2)]
    @test MOI.is_valid(optimizer, MOI.VariableIndex(2))
    @test !MOI.is_valid(optimizer, MOI.VariableIndex(3))

    MOI.optimize!(optimizer)

    # ensure that unmodified model is correct
    atol = config.atol
    rtol = config.rtol
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ -1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), v)) ≈ [1, 0] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[c]) ≈ -1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc1]) ≈ 0 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc2]) ≈ 1 atol=atol rtol=rtol

    # test default warm start
    itercold = optimizer.results.info.iter
    MOI.optimize!(optimizer)
    iterwarm = optimizer.results.info.iter
    @test iterwarm < itercold

    # test allocations
    allocs = @allocated MOI.optimize!(optimizer)
    @test allocs == 0

    # ensure that solving a second time results in the same answer after zeroing warm start
    zero_warm_start!(optimizer, values(idxmap.varmap), values(idxmap.conmap))
    test_optimizer_modification(m -> (), model, optimizer, idxmap, defaultoptimizer(), MOIT.TestConfig(atol=0.0, rtol=0.0))

    mapfrommodel(::MOI.AbstractOptimizer, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = idxmap[x]
    mapfrommodel(::MOI.ModelLike, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = x

    # change objective to min -2y
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        newobjf = MOI.ScalarAffineFunction([term(-2.0, mapfrommodel(m, y))], 0.0)
        F = typeof(newobjf)
        MOI.set(m, MOI.ObjectiveFunction{F}(), newobjf)
    end

    # add a constant to the objective
    objval_before = MOI.get(optimizer, MOI.ObjectiveValue())
    objconstant = 1.5
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        attr = MOI.ObjectiveFunction{Affine}()
        MOI.modify(m, attr, MOI.ScalarConstantChange(objconstant))
    end
    objval_after = MOI.get(optimizer, MOI.ObjectiveValue())
    @test objval_after ≈ objval_before + objconstant atol = 1e-8

    # change objective to min -y using ScalarCoefficientChange
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        attr = MOI.ObjectiveFunction{Affine}()
        MOI.modify(m, attr, MOI.ScalarCoefficientChange(mapfrommodel(m, y), -1.0))
    end
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ 0.5 * objval_before + objconstant atol = 1e-8

    # change x + y <= 1 to x + 2 y + 0.5 <= 1
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        attr = MOI.ConstraintFunction()
        ci =  mapfrommodel(m, c)
        MOI.set(m, attr, ci, MOI.ScalarAffineFunction(term.([1.0, 1.0, 1.0], mapfrommodel.(Ref(m), [x, x, y])), 0.5))
    end

    # change back to x + y <= 1 using ScalarCoefficientChange
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        ci = mapfrommodel(m, c)
        MOI.modify(m, ci, MOI.ScalarCoefficientChange(mapfrommodel(m, y), 1.0))
    end

    # flip the feasible set around from what it was originally and minimize +x
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        # objective
        newobjf = MOI.SingleVariable(mapfrommodel(m, x))
        F = typeof(newobjf)
        MOI.set(m, MOI.ObjectiveFunction{F}(), newobjf)

        # c
        attr = MOI.ConstraintFunction()
        ci = mapfrommodel(m, c)
        MOI.set(m, attr, ci, MOI.ScalarAffineFunction(term.([1.0, 1.0], mapfrommodel.(Ref(m), [x, y])), 0.0))

        attr = MOI.ConstraintSet()
        MOI.set(m, attr, ci, MOI.Interval(-1.0, Inf))

        # vc1
        attr = MOI.ConstraintSet()
        ci = mapfrommodel(m, vc1)
        MOI.set(m, attr, ci, MOI.Interval(-Inf, 0.))

        # vc2
        attr = MOI.ConstraintSet()
        ci = mapfrommodel(m, vc2)
        MOI.set(m, attr, ci, MOI.Interval(-Inf, 0.))
    end

    testflipped = function ()
        @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ -1 atol=atol rtol=rtol
        @test MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), v)) ≈ [-1, 0] atol=atol rtol=rtol
        @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[c]) ≈ 1 atol=atol rtol=rtol
        @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc1]) ≈ 0 atol=atol rtol=rtol
        @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc2]) ≈ -1 atol=atol rtol=rtol
    end
    testflipped()

    # update settings
    @test optimizer.results.info.status_polish == 0
    MOI.set(optimizer, OSQPSettings.Polish(), true)
    MOI.optimize!(optimizer)
    @test optimizer.results.info.status_polish == 1
    testflipped()
end

@testset "No CachingOptimizer: Vector problem modification after copy_to" begin
    # from basic.jl:
    model = OSQPModel{Float64}()
    x = MOI.add_variables(model, 2)
    P11 = 11.
    q = [3., 4.]
    u = [0., 0., -15, 100, 80]
    A = sparse(Float64[-1 0; 0 -1; -1 -3; 2 5; 3 4])
    I, J, coeffs = findnz(A)
    objf = MOI.ScalarQuadraticFunction(term.(q, x), [term(2 * P11, x[1], x[1])], 0.0)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objf)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    cf = MOI.VectorAffineFunction(MOI.VectorAffineTerm.(Int64.(I), term.(coeffs, map(j -> getindex(x, j), J))), -u)
    c = MOI.add_constraint(model, cf, MOI.Nonpositives(length(u)))

    optimizer = defaultoptimizer()
    idxmap = MOI.copy_to(optimizer, model)
    @test MOI.get(optimizer, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 2
    @test MOI.get(optimizer, MOI.ListOfVariableIndices()) == [MOI.VariableIndex(1), MOI.VariableIndex(2)]
    @test MOI.is_valid(optimizer, MOI.VariableIndex(2))
    @test !MOI.is_valid(optimizer, MOI.VariableIndex(3))

    MOI.optimize!(optimizer)

    # check result before modification
    atol = config.atol
    rtol = config.rtol
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ 20. atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), x)) ≈ [0.; 5.] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[c]) ≈ -[1.666666666666; 0.; 1.3333333; 0.; 0.] atol=atol rtol=rtol

    # test allocations
    allocs = @allocated MOI.optimize!(optimizer)
    @test allocs == 0

    mapfrommodel(::MOI.AbstractOptimizer, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = idxmap[x]
    mapfrommodel(::MOI.ModelLike, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = x

    # make random modifications to constraints
    randvectorconfig = MOIT.TestConfig(atol=Inf, rtol=1e-4)
    rng = MersenneTwister(1234)
    for i = 1 : 100
        newcoeffs = copy(coeffs)
        modindex = rand(rng, 1 : length(newcoeffs))
        newcoeffs[modindex] = 0
        newconst = 5 .* (rand(rng, length(u)) .- 0.5)
        test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), randvectorconfig) do m
            attr = MOI.ConstraintFunction()
            ci = mapfrommodel(m, c)
            newcf = MOI.VectorAffineFunction(MOI.VectorAffineTerm.(Int64.(I), term.(newcoeffs, map(j -> getindex(x, j), J))), newconst)
            MOI.set(m, attr, ci, newcf)
        end
    end
end

@testset "No CachingOptimizer: Warm starting" begin
    # define QP:
    # min 1/2 x'Px + q'x
    # s.t. Ax <= u   <=> -Ax + u in Nonnegatives
    #      Ax >= l   <=> Ax - l in Nonnegatives

    l = [1.; 0; 0];
    u = [1.; 0.7; 0.7]

    model = MOIU.UniversalFallback(OSQPModel{Float64}());
    optimizer = defaultoptimizer();

    x = MOI.add_variables(model, 2);
    objectiveFunction = MOI.ScalarQuadraticFunction{Float64}([MOI.ScalarAffineTerm(1.0, x[1]); MOI.ScalarAffineTerm(1.0, x[2])], [MOI.ScalarQuadraticTerm(4.0, x[1], x[1]); MOI.ScalarQuadraticTerm(1.0, x[1], x[2]); MOI.ScalarQuadraticTerm(2.0, x[2], x[2])] , 0);
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), objectiveFunction);
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE);
    Aneg = [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-1.0, x[1])),MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-1.0, x[2])),MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.0, x[1])),MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(-1.0, x[2]))];
    MOI.add_constraint(model, MOI.VectorAffineFunction(Aneg, u), MOI.Nonnegatives(3));
    A = [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[2])),MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, x[1])),MOI.VectorAffineTerm(3, MOI.ScalarAffineTerm(1.0, x[2]))];
    con1 = MOI.add_constraint(model, MOI.VectorAffineFunction(A, -u), MOI.Nonpositives(3));
    con2 = MOI.add_constraint(model, MOI.VectorAffineFunction(A, -l),MOI.Nonnegatives(3));


    MOI.empty!(optimizer);
    idxmap = MOI.copy_to(optimizer, model);
    MOI.optimize!(optimizer);

    # solve once to get optimal solution
    x_sol = MOI.get(optimizer, MOI.VariablePrimal(), getindex.(Ref(idxmap), x));
    y_c1 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con1]);
    y_c2 = MOI.get(optimizer, MOI.ConstraintDual(), idxmap[con2]);
    y_c1_rows = OSQP.MathOptInterfaceOSQP.constraint_rows(optimizer.rowranges, idxmap[con1]);
    y_c2_rows = OSQP.MathOptInterfaceOSQP.constraint_rows(optimizer.rowranges, idxmap[con2]);

    # provide warm start values to the model
    MOI.set.(model, MOI.VariablePrimalStart(), x, x_sol)
    MOI.set.(model, MOI.ConstraintDualStart(), [con1, con2], [y_c1, y_c2])
    MOI.empty!(optimizer);
    idxmap = MOI.copy_to(optimizer, model);

    # check that internal variables are set correctly
    @test optimizer.warmstartcache.x.data == x_sol
    @test optimizer.warmstartcache.y.data[y_c1_rows] == -y_c1
    @test optimizer.warmstartcache.y.data[y_c2_rows] == -y_c2
end

@testset "Vector equality constraint" begin
    # Minimize ||A x - b||^2 = x' A' A x - (2 * A' * b)' x + b' * b
    # subject to C x = d

    generate_problem_data = function (rng, n, m)
        A = rand(rng, n, n)
        b = rand(rng, n)
        C = rand(rng, m, n)
        d = rand(rng, m)
        C⁺ = pinv(C)
        Q = I - C⁺ * C
        expected = Q * (pinv(A * Q) * (b - A * C⁺ * d)) + C⁺ * d # note: can be quite badly conditioned
        @test C * expected ≈ d atol = 1e-12

        P = Symmetric(sparse(triu(A' * A)))
        q = -2 * A' * b
        r = b' * b

        A, b, C, d, P, q, r, expected
    end

    make_objective = function (P, q, r, x)
        I, J, coeffs = findnz(P.data)
        MOI.ScalarQuadraticFunction(term.(q, x), term.(2 * coeffs, map(i -> x[i]::MOI.VariableIndex, I), map(j -> x[j]::MOI.VariableIndex, J)), r)
    end

    make_constraint_fun = function (C, d, x)
        I, J, coeffs = findnz(sparse(C))
        cf = MOI.VectorAffineFunction(MOI.VectorAffineTerm.(Int64.(I), term.(coeffs, map(j -> getindex(x, j)::MOI.VariableIndex, J))), -d)
    end

    check_results = function (optimizer, idxmap, x, A, b, expected)
        @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get.(Ref(optimizer), Ref(MOI.VariablePrimal()), getindex.(Ref(idxmap), x)) ≈ expected atol = 1e-4
        @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ norm(A * expected - b)^2 atol = 1e-4
    end

    n = 8
    m = 2
    rng = MersenneTwister(1234)

    A, b, C, d, P, q, r, expected = generate_problem_data(rng, n, m)
    model = OSQPModel{Float64}()
    x = MOI.add_variables(model, n)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), make_objective(P, q, r, x))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    c = MOI.add_constraint(model, make_constraint_fun(C, d, x), MOI.Zeros(length(d)))

    optimizer = defaultoptimizer()
    idxmap = MOI.copy_to(optimizer, model)
    MOI.optimize!(optimizer)
    check_results(optimizer, idxmap, x, A, b, expected)

    x = [idxmap[xi] for xi in x]
    for i = 1 : 10
        A, b, C, d, P, q, r, expected = generate_problem_data(rng, n, m)
        MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), make_objective(P, q, r, x))
        attr = MOI.ConstraintFunction()
        MOI.set(optimizer, attr, idxmap[c], make_constraint_fun(C, d, x))
        attr = MOI.ConstraintSet()
        MOI.set(optimizer, attr, idxmap[c], MOI.Zeros(length(d))) # noop, but ok
        MOI.optimize!(optimizer)
        check_results(optimizer, idxmap, x, A, b, expected)
    end
end

@testset "RawSolver" begin
    optimizer = defaultoptimizer()
    let inner = MOI.get(optimizer, MOI.RawSolver())
        @test inner.workspace == C_NULL
    end

    @test MOI.get(optimizer, MOI.SolverName()) == "OSQP"

    model = OSQPModel{Float64}()
    MOI.empty!(model)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, x, MOI.GreaterThan(2.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))

    MOI.copy_to(optimizer, model)
    inner = MOI.get(optimizer, MOI.RawSolver())
    @test inner.workspace != C_NULL
end

# TODO: consider moving to MOIT. However, current default copy_to is fine with BadObjectiveModel.
struct BadObjectiveModel <: MOIT.BadModel end # objective sense is not FEASIBILITY_SENSE, but can't get objective function
MOI.get(src::BadObjectiveModel, ::MOI.ObjectiveSense) = MOI.MIN_SENSE
struct ExoticFunction <: MOI.AbstractScalarFunction end
MOI.get(src::BadObjectiveModel, ::MOI.ObjectiveFunctionType) = ExoticFunction

@testset "failcopy" begin
    # TODO change OSQPOptimizer() to OSQP.Optimizer() in OSQP v0.6
    optimizer = OSQPOptimizer()
    MOIT.failcopytestc(optimizer)
    @test_throws MOI.UnsupportedAttribute{MOI.ObjectiveFunction{ExoticFunction}} MOI.copy_to(optimizer, BadObjectiveModel())
end
