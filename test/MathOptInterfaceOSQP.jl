using OSQP.MathOptInterfaceOSQP
using Base.Test

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

using MathOptInterfaceUtilities
const MOIU = MathOptInterfaceUtilities

MOIU.@instance(OSQPInstanceData, # instancename
    (), # scalarsets
    (Interval,), # typedscalarsets
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
    deleteat!(excludes, findfirst(excludes, "linear10"))

    MOIT.contlineartest(MOIU.InstanceManager(OSQPInstanceData{Float64}(), OSQPInstance()), config, excludes)
end
