__precompile__()

module OSQP

export OSQPMathProgBaseInterface

# Compatibility stuff
using Compat
using Compat.SparseArrays
using Compat.Iterators




if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("OSQP not properly installed. Please run Pkg.build(\"OSQP\")")
end


function __init__()
    # Get version
    ver_array = split(version(), ".")
    ver_string = string(ver_array[1], ".", ver_array[2], ".", ver_array[3])  # Get string without dev vers
    vnum = VersionNumber(ver_string)


    depsdir = realpath(joinpath(dirname(@__FILE__), "..", "deps"))
    if (vnum.major != 0 && vnum.minor != 2)
        error("Current OSQP version installed is $(osqp_version()), but we require version 0.2.*. Delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"OSQP\").")
    end
end


include("constants.jl")
include("types.jl")
include("interface.jl")
include("mpbinterface.jl")

end # module
