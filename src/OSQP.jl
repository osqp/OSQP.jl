__precompile__()

module OSQP

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("OSQP not properly installed. Please run Pkg.build(\"OSQP\")")
end



import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen



function __init__()
    vnum = VersionNumber(version())
    # default binaries need access to Julia's lapack symbols
    if is_unix()
        dlopen(Base.liblapack_name, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
    end

    # TODO: Add version check!
    # depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
    # if !(vnum.major == 1 && vnum.minor in [1,2])
    #     error("Current OSQP version installed is $(osqp_version()), but we require version 1.1.* or 1.2.*. On Linux and Windows, delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"OSQP\"). On OS X, run `using Homebrew; Homebrew.update()` in Julia.")
    # end
    # if is_linux() && vnum.major == 1 && vnum.minor == 1
    #     warn("Current OSQP version installed is $(osqp_version()), but a newer version (1.2.*) is available. To upgrade, delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"OSQP\").")
    # end
end





include("types.jl")
include("interface.jl")
include("utils.jl")

end # module
