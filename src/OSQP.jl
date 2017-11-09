__precompile__()

module OSQP

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("OSQP not properly installed. Please run Pkg.build(\"OSQP\")")
end



import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen



function __init__()
    # Get version number
    # vnum = VersionNumber(version())

    # Get version with dev inside (debug)
    vnum = VersionNumber(version()[1:end-5])

    depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
    if (vnum.major != 0 && vnum.minor != 2)
	error("Current OSQP version installed is $(osqp_version()), but we require version 0.2.*. Delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"OSQP\").")
    end
end





include("constants.jl")
include("types.jl")
include("interface.jl")

end # module
