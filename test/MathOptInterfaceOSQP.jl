using OSQP.MathOptInterfaceOSQP
using Base.Test

using MathOptInterface
const MOI = MathOptInterface

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities

MOIU.@model(OSQPCachingOptimizer, # modelname
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
    MOIT.contlineartest(MOIU.CachingOptimizer(OSQPCachingOptimizer{Float64}(), optimizer), config, excludes)
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
    MOIT.contquadratictest(MOIU.CachingOptimizer(OSQPCachingOptimizer{Float64}(), optimizer), config, excludes)
end
