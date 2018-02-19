using Compat # for CartesianIndices
using OSQP.MathOptInterfaceOSQP
using Base.Test

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterface.Test
const MOIT = MathOptInterface.Test

using MathOptInterface.Utilities
const MOIU = MathOptInterface.Utilities

@testset "ProblemModificationCache" begin
    rng = MersenneTwister(1234)
    n = 15
    m = 10
    q = randn(rng, n)
    P = (X = sprandn(rng, n, n, 0.1); X' * X)
    P += eps() * Base.I # needed for a test later on
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
    @test qmod_update_results.x ≈ qmod_setup_results.x atol = 1e-8

    # Modify A, ensure that updating results in the same solution as calling setup! with the modified A and q
    (I, J) = findn(A)
    Amodindex = rand(rng, 1 : nnz(A))
    row = I[Amodindex]
    col = J[Amodindex]
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
    @test Amod_update_results.x ≈ Amod_setup_results.x atol = 1e-8

    # MatrixModificationCache: colon indexing
    modcache.P[:] = 0.
    for i = 1 : n
        modcache.P[i, i] = 1.
    end
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    Pmod_update_results = OSQP.solve!(model)
    model4 = OSQP.Model()
    Pmod = speye(n, n)
    OSQP.setup!(model4; P = Pmod, q = modcache.q.data, A = Amod, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    Pmod_setup_results = OSQP.solve!(model4)
    @test Pmod_update_results.x ≈ Pmod_setup_results.x atol = 1e-8

    # Modifying the sparsity pattern is not allowed
    nzinds = map(CartesianIndex, zip(I, J))
    zinds = zinds = setdiff(vec(CartesianIndices(A)), nzinds)
    for I in zinds
        @test_throws ArgumentError modcache.A[I[1], I[2]] = randn(rng)
    end
    @test_throws ArgumentError modcache.A[:] = 1
end

# FIXME: type piracy. Generalize and move to MOIU.
MOI.canget(optimizer::MOIU.CachingOptimizer, ::MOI.ConstraintPrimal, ::Type{<:MOI.ConstraintIndex}) = true
function MOI.get(optimizer::MOIU.CachingOptimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{<:MOI.SingleVariable, <:Any})
    f = MOI.get(optimizer, MOI.ConstraintFunction(), ci)
    MOI.get(optimizer, MOI.VariablePrimal(), f.variable)
end
function MOI.get(optimizer::MOIU.CachingOptimizer, ::MOI.ConstraintPrimal, ci::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction, <:Any})
    f = MOI.get(optimizer, MOI.ConstraintFunction(), ci)
    ret = f.constant
    n = length(f.variables)
    length(f.coefficients) == n || error()
    for i = 1 : n
        ret += f.coefficients[i] * MOI.get(optimizer, MOI.VariablePrimal(), f.variables[i])
    end
    ret
end

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

function defaultoptimizer()
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOI.set!(optimizer, OSQPSettings.MaxIter(), 10000)
    MOI.set!(optimizer, OSQPSettings.CheckTermination(), true)
    optimizer
end

@testset "CachingOptimizer: linear problems" begin
    excludes = String[]
    push!(excludes, "linear7") # vector constraints
    optimizer = defaultoptimizer()
    MOIT.contlineartest(MOIU.CachingOptimizer(OSQPModel{Float64}(), optimizer), config, excludes)
end

@testset "CachingOptimizer: quadratic problems" begin
    excludes = String[]
    push!(excludes, "quadratic4") # QCP
    push!(excludes, "quadratic5") # QCP
    push!(excludes, "quadratic6") # QCP
    push!(excludes, "quadratic7") # SOCP
    optimizer = defaultoptimizer()
    MOIT.contquadratictest(MOIU.CachingOptimizer(OSQPModel{Float64}(), optimizer), config, excludes)
end

function test_optimizer_modification(modfun::Base.Callable, model::MOI.ModelLike, optimizer::T, idxmap::MOIU.IndexMap,
        cleanoptimizer::T, config::MOIT.TestConfig) where T<:MOI.AbstractOptimizer
    # apply modfun to both the model and the optimizer
    modfun(model)
    modfun(optimizer)

    # copy model into clean optimizer
    copyresult = MOI.copy!(cleanoptimizer, model)
    @test copyresult.status == MOI.CopySuccess
    @test copyresult.indexmap.varmap == idxmap.varmap
    @test copyresult.indexmap.conmap == idxmap.conmap

    # call optimize! on both optimizers
    MOI.optimize!(optimizer)
    MOI.optimize!(cleanoptimizer)

    # compare results
    atol = config.atol
    rtol = config.rtol
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.get(cleanoptimizer, MOI.TerminationStatus())
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.get(cleanoptimizer, MOI.PrimalStatus())
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ MOI.get(cleanoptimizer, MOI.ObjectiveValue()) atol=atol rtol=rtol

    modelvars = MOI.get(model, MOI.ListOfVariableIndices())
    for v_model in modelvars
        v_optimizer = idxmap[v_model]
        @test MOI.get(optimizer, MOI.VariablePrimal(), v_optimizer) ≈ MOI.get(cleanoptimizer, MOI.VariablePrimal(), v_optimizer) atol=atol rtol=rtol
    end

    if config.duals
        @test MOI.get(optimizer, MOI.DualStatus()) == MOI.get(cleanoptimizer, MOI.DualStatus())
        for (F, S) in MOI.get(model, MOI.ListOfConstraints())
            cis_model = MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
            for ci_model in cis_model
                ci_optimizer = idxmap[ci_model]
                @test MOI.get(optimizer, MOI.ConstraintDual(), ci_optimizer) ≈ MOI.get(cleanoptimizer, MOI.ConstraintDual(), ci_optimizer) atol=atol rtol=rtol
            end
        end
    end
end

function zero_warm_start!(optimizer::MOI.ModelLike, vars, cons)
    for vi in vars
        MOI.set!(optimizer, MOI.VariablePrimalStart(), vi, 0.0)
    end
    for ci in cons
        MOI.set!(optimizer, MOI.ConstraintDualStart(), ci, 0.0)
    end
end

@testset "No CachingOptimizer: problem modification after copy!" begin
    # Initial setup: modified version of MOIT.linear1test
    # min -x
    # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
    #       x, y >= 0   (x, y ∈ Nonnegatives)
    model = OSQPModel{Float64}()
    MOI.empty!(model)
    v = MOI.addvariables!(model, 2)
    x, y = v
    cf = MOI.ScalarAffineFunction([v; v; v], [0.0,0.0,1.0,1.0,0.0,0.0], 0.0)
    c = MOI.addconstraint!(model, cf, MOI.Interval(-Inf, 1.0))
    vc1 = MOI.addconstraint!(model, MOI.SingleVariable(v[1]), MOI.Interval(0.0, Inf))
    vc2 = MOI.addconstraint!(model, v[2], MOI.Interval(0.0, Inf))
    objf = MOI.ScalarAffineFunction([v; v; v], [0.0,0.0,-1.0,0.0,0.0,0.0], 0.0)
    MOI.set!(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), objf)
    MOI.set!(model, MOI.ObjectiveSense(), MOI.MinSense)

    optimizer = defaultoptimizer()
    copyresult = MOI.copy!(optimizer, model)
    idxmap = copyresult.indexmap

    # Since OSQP automatically warm starts, we need to manually zero to get repeatable results
    zero_warm_start!(optimizer, values(idxmap.varmap), values(idxmap.conmap))
    MOI.optimize!(optimizer)

    # ensure that unmodified model is correct
    atol = config.atol
    rtol = config.rtol
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ -1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), getindex.(idxmap, v)) ≈ [1, 0] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[c]) ≈ -1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc1]) ≈ 0 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc2]) ≈ 1 atol=atol rtol=rtol

    # ensure that solving a second time results in the same answer (after zeroing warm start)
    zero_warm_start!(optimizer, values(idxmap.varmap), values(idxmap.conmap))
    # TODO: I would have expected this to pass:
    # test_optimizer_modification(m -> (), model, optimizer, idxmap, defaultoptimizer(), MOIT.TestConfig(atol=0.0, rtol=0.0))
    test_optimizer_modification(m -> (), model, optimizer, idxmap, defaultoptimizer(), config)

    mapfrommodel(::MOI.AbstractOptimizer, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = idxmap[x]
    mapfrommodel(::MOI.ModelLike, x::Union{MOI.VariableIndex, <:MOI.ConstraintIndex}) = x

    # change objective to min -2y
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        newobjf = MOI.ScalarAffineFunction([mapfrommodel(m, y)], [-2.0], 0.0)
        F = typeof(newobjf)
        @test MOI.canset(m, MOI.ObjectiveFunction{F}())
        MOI.set!(m, MOI.ObjectiveFunction{F}(), newobjf)
    end

    # change x + y <= 1 to x + 2 y <= 1
    zero_warm_start!(optimizer, values(idxmap.varmap), values(idxmap.conmap))
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        @test MOI.canmodifyconstraint(m, mapfrommodel(m, c), MOI.ScalarAffineFunction{Float64})
        MOI.modifyconstraint!(m, mapfrommodel(m, c), MOI.ScalarAffineFunction(mapfrommodel.(m, [x, x, y]), [1.0, 1.0, 1.0], 0.0))
    end

    # flip the feasible set around from what it was originally and minimize +x
    test_optimizer_modification(model, optimizer, idxmap, defaultoptimizer(), config) do m
        # objective
        newobjf = MOI.SingleVariable(mapfrommodel(m, x))
        F = typeof(newobjf)
        @test MOI.canset(m, MOI.ObjectiveFunction{F}())
        MOI.set!(m, MOI.ObjectiveFunction{F}(), newobjf)

        # c
        @test MOI.canmodifyconstraint(m, mapfrommodel(m, c), MOI.ScalarAffineFunction{Float64})
        MOI.modifyconstraint!(m, mapfrommodel(m, c), MOI.ScalarAffineFunction(mapfrommodel.(m, [x, y]), [1.0, 1.0], 0.0))
        @test MOI.canmodifyconstraint(m, mapfrommodel(m, c), MOI.Interval{Float64})
        MOI.modifyconstraint!(m, mapfrommodel(m, c), MOI.Interval(-1.0, Inf))

        # vc1
        @test MOI.canmodifyconstraint(m, mapfrommodel(m, vc1), MOI.Interval{Float64})
        MOI.modifyconstraint!(m, mapfrommodel(m, vc1), MOI.Interval(-Inf, 0.))

        # vc2
        @test MOI.canmodifyconstraint(m, mapfrommodel(m, vc2), MOI.Interval{Float64})
        MOI.modifyconstraint!(m, mapfrommodel(m, vc2), MOI.Interval(-Inf, 0.))
    end

    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ -1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.VariablePrimal(), getindex.(idxmap, v)) ≈ [-1, 0] atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[c]) ≈ 1 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc1]) ≈ 0 atol=atol rtol=rtol
    @test MOI.get(optimizer, MOI.ConstraintDual(), idxmap[vc2]) ≈ -1 atol=atol rtol=rtol
end

# TODO: consider moving to MOIT. However, current defaultcopy! is fine with BadObjectiveModel.
struct BadObjectiveModel <: MOIT.BadModel end # objective sense is not FeasibilitySense, but can't get objective function
MOI.canget(src::BadObjectiveModel, ::MOI.ObjectiveSense) = true
MOI.get(src::BadObjectiveModel, ::MOI.ObjectiveSense) = MOI.MinSense
MOI.canget(src::BadObjectiveModel, ::MOI.ObjectiveFunction{<:Any}) = false

@testset "failcopy" begin
    optimizer = OSQPOptimizer()
    MOIT.failcopytestc(optimizer)
    MOIT.failcopytest(optimizer, BadObjectiveModel(), MOI.CopyOtherError)
end
