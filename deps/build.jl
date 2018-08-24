using BinDeps
using Compat

# Libdl has been moved to a standard library module
using Compat.Libdl

@BinDeps.setup

# Add library dependency for direct method
osqp = library_dependency("osqp", aliases=["libosqp"])

# Current version
version = "0.4.0"
qdldl_version = "0.1.2"

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
provides(Sources, URI("https://github.com/oxfordcontrol/osqp/archive/v$version.tar.gz"),
    [osqp], unpacked_dir="osqp-$version", os = :Linux)

# Define directories locations deps/usr and deps/src
prefix = joinpath(BinDeps.depsdir(osqp),"usr")
srcdir = joinpath(BinDeps.depsdir(osqp),"src","osqp-$version")
blddir = joinpath(srcdir, "build")
downloadsdir = joinpath(BinDeps.depsdir(osqp),"downloads")

# Define library name
libname = "libosqp.$(dlext)"

# File move rule to move QDLDL sources
struct FileMoveRule <: BinDeps.BuildStep
    src::AbstractString
    dest::AbstractString
end
Base.run(fc::FileMoveRule) = isfile(fc.dest) || mv(fc.src, fc.dest, force=true)

provides(SimpleBuild,
    (@build_steps begin
    GetSources(osqp)
    # Copy QDLDL sources inside OSQP sources
    FileDownloader("https://github.com/oxfordcontrol/qdldl/archive/v$qdldl_version.tar.gz",
		   joinpath(downloadsdir, "external", "qdldl-$qdldl_version.tar.gz"))
    FileUnpacker(joinpath(downloadsdir, "external", "qdldl-$qdldl_version.tar.gz"),
		 joinpath(srcdir, downloadsdir),
		 "qdldl-$qdldl_version")
    FileMoveRule(joinpath(srcdir, downloadsdir, "qdldl-$qdldl_version"),
		 joinpath(srcdir, "lin_sys", "direct", "qdldl", "qdldl_sources"))
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
