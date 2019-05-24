using BinaryProvider

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

# Current version
version = "0.5.0"

# Get current operating system
osqp_platform =
if Sys.islinux()
    "linux"
elseif Sys.isapple()
    "mac"
elseif Sys.iswindows()
    "windows"
else
    error("Platform not supported!")
end

subfolder="osqp-$version-$osqp_platform$(Sys.WORD_SIZE)"
products = [
        LibraryProduct(joinpath(prefix, subfolder, "lib"), ["libosqp"], :osqp),
    ]

# Provide binaries for each operating system

function get_hash_code(platform)
    archive_name = "osqp-$version-$platform.tar.gz"
    readlines(download("https://dl.bintray.com/bstellato/generic/OSQP/$version/$archive_name.sha256"))[1]
end

function get_url(platform)
    archive_name = "osqp-$version-$platform"
    return "https://dl.bintray.com/bstellato/generic/OSQP/$version/$archive_name.tar.gz"
end

download_info = Dict(
    Linux(:x86_64) => (get_url("linux64"), get_hash_code("linux64")),
    MacOS(:x86_64) => (get_url("mac64"), get_hash_code("mac64")),
    Windows(:x86_64) => (get_url("windows64"), get_hash_code("windows64")),
    Windows(:i686) => (get_url("windows32"), get_hash_code("windows32")),
)

# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)
if haskey(download_info, platform_key())
    url, tarball_hash = download_info[platform_key()]

    # Check if this build.jl is providing new versions of the binaries, and
    # if so, ovewrite the current binaries even if they were installed by the user
    if unsatisfied || !isinstalled(url, tarball_hash; prefix=prefix)
        # Download and install binaries
        install(url, tarball_hash; prefix=prefix, force=true, verbose=verbose)
    end
elseif unsatisfied
    # If we don't have a BinaryProvider-compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform $(triplet(platform_key_abi())) is not supported by this package!")
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
