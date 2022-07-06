"""
    OSQPBackend

The supertype for all backend library implementations of OSQP to implement
a library pointer type underneath.
"""
abstract type OSQPBackend end


"""
    osqp_ccall(func, b::backend, returnType, paramTypes, args...) where {backend <: OSQPBackend}

Call a function `func` inside the OSQP backend library specified by `backend`. The remaining
arguments are the same as a normal ccall.
"""
@generated function osqp_ccall(func, b::backend, returnType, paramTypes, args...) where {backend <: OSQPBackend}
    return _osqp_ccall( func, b, returnType, paramTypes, args...)
end

function _osqp_ccall( func, ::Type{backend}, returnType, paramTypes, args...) where {backend <: OSQPBackend}
    libname = getlibrary(backend())

    return :(ccall((func, $(libname)), returnType, paramTypes, args...))
end

export OSQPBuiltin


"""
    OSQPBuiltin

Use the builtin OSQP linear algebra functions and linear system solvers.
"""
struct OSQPBuiltin <: OSQPBackend end

function getlibrary( ::OSQPBuiltin )
    return osqp_builtin
end
