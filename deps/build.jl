using BinDeps

@BinDeps.setup

# Add library dependency for direct method
osqp = library_dependency("osqp", aliases=["libosqp"])

# Current version
version = "0.2.0.dev5"

if Sys.KERNEL == :Linux
	sys_name_str = "linux"
elseif Sys.KERNEL == :Darwin
	sys_name_str = "mac"
else
	sys_name_str = "windows"
end


# Provide binaries for each operative system
archive_name="osqp-$version-$(sys_name_str)$(Sys.WORD_SIZE)"

provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/$archive_name.tar.gz"), [osqp], unpacked_dir="$archive_name/lib")

# # Windows
# provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/osqp-$version-windows$(Sys.WORD_SIZE).tar.gz"), [osqp], os = :Windows, unpacked_dir="osqp-$version")
#
# # Mac
# provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/osqp-$version-mac64.tar.gz"), [osqp], os = :Darwin, unpacked_dir="osqp-$version")
#
#
# # Linux
# provides(Binaries, URI("https://dl.bintray.com/bstellato/generic/OSQP/$version/osqp-$version-linux64.tar.gz"), [osqp], os = :Linux, unpacked_dir="osqp-$version")


# https://dl.bintray.com/bstellato/generic/OSQP/0.2.0.dev5/:osqp-0.2.0.dev5-windows32.tar.gz
# provides(Sources, URI("https://github.com/oxfordcontrol/osqp/archive/$version.tar.gz"),
    # [osqp], unpacked_dir="osqp-$version")

# provides(Sources, URI("https://github.com/oxfordcontrol/osqp/archive/v$version.tar.gz"),
#     [osqp], unpacked_dir="osqp-$version")

# Define directories locations deps/usr and deps/src
# prefix = joinpath(BinDeps.depsdir(osqp),"usr")
# srcdir = joinpath(BinDeps.depsdir(osqp),"src","osqp-$version")
# blddir = joinpath(srcdir, "build")

# # Define library name
# libname = "libosqp.$(Libdl.dlext)"
#
#
# provides(SimpleBuild,
#     (@build_steps begin
#         GetSources(osqp)
#         CreateDirectory(joinpath(prefix, "lib"))
#         @build_steps begin
#             # ChangeDirectory(srcdir)
#             CreateDirectory(blddir)
#             @build_steps begin
#                 ChangeDirectory(blddir)
#                 FileRule(joinpath(prefix, "lib", libname), @build_steps begin
#                     `cmake -G "Unix Makefiles" -DUNITTESTS=OFF ..`
#                     `make osqp`
#                     `mv out/$libname $prefix/lib`
#                 end)
#             end
#         end
#     end),
# [osqp])
#
#
#
#
# TODO: Add proper dependencies download
# if is_apple()
#     using Homebrew
#     provides(Homebrew.HB, "osqp", osqp, os = :Darwin)
# end


# TODO: Provide sources only for Unix
# provides(Sources, URI("https://github.com/oxfordcontrol/osqp/archive/v$version.tar.gz"),
#     [osqp], os = :Unix, unpacked_dir="osqp-$version")




@BinDeps.install Dict(:osqp => :osqp)
