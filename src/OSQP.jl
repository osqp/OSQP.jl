module OSQP

export OSQPMathProgBaseInterface
using SparseArrays
using LinearAlgebra

if VERSION < v"1.3"
    if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        include("../deps/deps.jl")
    else
        error("OSQP not properly installed. Please run Pkg.build(\"OSQP\")")
    end
else
    using OSQP_jll
end


function __init__()
    # Get version
    if VERSION < v"1.3"
        ver_array = split(version(), ".")
        ver_string = string(ver_array[1], ".", ver_array[2], ".", ver_array[3])  # Get string without dev vers
        vnum = VersionNumber(ver_string)


        depsdir = realpath(joinpath(dirname(@__FILE__), "..", "deps"))
        if (vnum.major != 0 && vnum.minor != 6)
            error("Current OSQP version installed is $(osqp_version()), but we require version 0.6.*. Delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"OSQP\").")
        end
    end
end

include("constants.jl")
include("types.jl")
include("interface.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")
const Optimizer = MathOptInterfaceOSQP.Optimizer

end # module
