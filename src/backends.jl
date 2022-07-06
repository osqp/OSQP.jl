"""
    OSQPBackend


"""
abstract type OSQPBackend end

@generated function osqp_ccall(func, b::backend, args...) where {backend <: OSQPBackend}
    libname = getlibrary(backend())

    return :(ccall((func, $(libname)), args...))
    #return _osqp_ccall(func, b, args...)
end

function _osqp_ccall( func, ::Type{backend}, args...) where {backend <: OSQPBackend}
    libname = getlibrary(backend())

    return :(ccall((func, $(libname)), args...))
end

"""
    OSQPBuiltin

Use the builtin OSQP linear algebra functions and linear system solvers.
"""
struct OSQPBuiltin <: OSQPBackend end

function getlibrary( ::OSQPBuiltin )
    return osqp_builtin
end
