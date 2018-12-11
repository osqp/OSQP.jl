using BinDeps
using Compat

# Libdl has been moved to a standard library module
using Compat.Libdl

@BinDeps.setup

# Add library dependency for direct method
osqp = library_dependency("osqp", aliases=["libosqp"])

# Current version
version = "0.5.0"

# Get current operating system
osqp_platform =
if Compat.Sys.islinux()
    "linux"
elseif Compat.Sys.isapple()
    "mac"
elseif Compat.Sys.iswindows()
    "windows"
else
    error("Platform not supported!")
end

# Provide binaries for each operating system
archive_name="osqp-$version-$osqp_platform$(Sys.WORD_SIZE)"
provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/$archive_name.tar.gz"), [osqp], unpacked_dir="$archive_name/lib", os=:Darwin)
provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/$archive_name.tar.gz"), [osqp], unpacked_dir="$archive_name/lib", os=:Windows)

# Build from sources on Linux
provides(Sources, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/osqp-$version.tar.gz"), [osqp], unpacked_dir="osqp-$version", os = :Linux)

# Define directories locations deps/usr and deps/src
prefix = joinpath(BinDeps.depsdir(osqp),"usr")
srcdir = joinpath(BinDeps.depsdir(osqp),"src","osqp-$version")
blddir = joinpath(srcdir, "build")

# Define library name
libname = "libosqp.$(dlext)"

provides(SimpleBuild,
    (@build_steps begin
    GetSources(osqp)
    CreateDirectory(joinpath(prefix, "lib"))
    @build_steps begin
        CreateDirectory(blddir)
        @build_steps begin
        ChangeDirectory(blddir)
        FileRule(joinpath(prefix, "lib", libname), @build_steps begin
            `cmake -G "Unix Makefiles" ..`
            `make osqp`
            `mv out/$libname $prefix/lib`
        end)
        end
    end
    end),
[osqp], os = :Linux)


@BinDeps.install Dict(:osqp => :osqp)
