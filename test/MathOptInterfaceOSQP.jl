using OSQP.MathOptInterfaceOSQP
using Base.Test

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities

@testset "ProblemModificationCache" begin
    rng = MersenneTwister(1234)
    n = 15
    m = 10
    q = randn(rng, n)
    P = (X = sprandn(rng, n, n, 0.1); X' * X)
    l = -rand(rng, m)
    u = rand(rng, m)
    A = sprandn(rng, m, n, 0.6)
    modcache = MathOptInterfaceOSQP.ProblemModificationCache(P, q, A, l, u)

    model = OSQP.Model()
    OSQP.setup!(model; P = P, q = q, A = A, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    baseresults = OSQP.solve!(model)
    @test baseresults.info.status == :Solved

    # VectorModificationCache basics
    @test !modcache.qcache.dirty
    modcache.qcache[3] = 5.
    @test modcache.qcache.dirty
    @test modcache.qcache.data[3] == 5.

    # Process q modifications, ensure that updating results in the same solution as calling setup! with the modified q
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    @test !modcache.qcache.dirty
    qmod_update_results = OSQP.solve!(model)
    @test !isapprox(baseresults.x, qmod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    model2 = OSQP.Model()
    OSQP.setup!(model2; P = P, q = modcache.qcache.data, A = A, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    qmod_setup_results = OSQP.solve!(model2)
    @test qmod_update_results.x ≈ qmod_setup_results.x atol = 1e-8

    # MatrixModificationCache basics
    (I, J) = findn(A)
    Amodindex = rand(rng, 1 : nnz(A))
    row = I[Amodindex]
    col = J[Amodindex]
    val = randn(rng)
    modcache.Acache[row, col] = val
    @test any(x -> x == val, modcache.Acache.modifications)

    # Process A modifications, ensure that updating results in the same solution as calling setup! with the modified A and q
    MathOptInterfaceOSQP.processupdates!(model, modcache)
    @test all(iszero, modcache.Acache.modifications)
    Amod_update_results = OSQP.solve!(model)
    @test !isapprox(baseresults.x, Amod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    @test !isapprox(qmod_update_results.x, Amod_update_results.x; atol = 1e-1) # ensure that new results are significantly different
    Amod = copy(A)
    Amod[row, col] = val
    model3 = OSQP.Model()
    OSQP.setup!(model3; P = P, q = modcache.qcache.data, A = Amod, l = l, u = u, verbose = false, eps_abs = 1e-8, eps_rel = 1e-16)
    Amod_setup_results = OSQP.solve!(model3)
    @test Amod_update_results.x ≈ Amod_setup_results.x atol = 1e-8

    # Modifying the sparsity pattern is not allowed
    nzinds = map(CartesianIndex, zip(I, J))
    zinds = zinds = setdiff(vec(CartesianIndices(A)), nzinds)
    for I in zinds
        @test_throws ArgumentError modcache.Acache[I[1], I[2]] = randn(rng)
    end
end

MOIU.@model(OSQPCacheModel, # modelname
    (), # scalarsets
    (Interval, LessThan, GreaterThan, EqualTo), # typedscalarsets
    (), # vectorsets
    (), # typedvectorsets
    (SingleVariable,), # scalarfunctions
    (ScalarAffineFunction, ScalarQuadraticFunction), # typedscalarfunctions
    (), # vectorfunctions
    () # typedvectorfunctions
)

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

@testset "Continuous linear problems" begin
    excludes = String[]
    push!(excludes, "linear7") # vector constraints
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOIT.contlineartest(MOIU.CachingOptimizer(OSQPCacheModel{Float64}(), optimizer), config, excludes)
end

@testset "Continuous quadratic problems" begin
    excludes = String[]
    push!(excludes, "quadratic4") # QCP
    push!(excludes, "quadratic5") # QCP
    push!(excludes, "quadratic6") # QCP
    push!(excludes, "quadratic7") # SOCP
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-16)
    MOIT.contquadratictest(MOIU.CachingOptimizer(OSQPCacheModel{Float64}(), optimizer), config, excludes)
end
