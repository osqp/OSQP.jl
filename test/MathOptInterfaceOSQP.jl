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

const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Continuous linear problems" begin
    excludes = collect(keys(MOIT.contlineartests))
    deleteat!(excludes, findfirst(excludes, "linear1"))
    optimizer = OSQPOptimizer()
    MOI.set!(optimizer, OSQPSettings.Verbose(), false)
    MOI.set!(optimizer, OSQPSettings.EpsAbs(), 1e-8)
    MOI.set!(optimizer, OSQPSettings.EpsRel(), 1e-8)
    MOIT.contlineartest(MOIU.CachingOptimizer(OSQPCachingOptimizer{Float64}(), optimizer), config, excludes)
end
