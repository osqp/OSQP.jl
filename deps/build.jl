using BinaryProvider

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

# Current version
version = "0.6.0"

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
dl_info = choose_download(download_info, platform_key_abi())
if dl_info === nothing && unsatisfied
    # If we don't have a compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
end

# If we have a download, and we are unsatisfied (or the version we're
# trying to install is not itself installed) then load it up!
if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
    # Download and install binaries
    install(dl_info...; prefix=prefix, force=true, verbose=verbose)
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
