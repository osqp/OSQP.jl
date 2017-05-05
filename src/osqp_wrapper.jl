# Wrapper for the low level functions defined in https://github.com/oxfordcontrol/osqp/blob/master/include/osqp.h

# macro to call a OSQP C function
macro osqp_ccall(func, args...)
    f = "osqp_$(func)"
    args = map(esc,args)
    is_unix() && return quote
        ccall(($f,OSQP.osqp), $(args...))
    end
    is_windows() && VERSION < v"0.6-" && return quote
        ccall(($f,OSQP.osqp), stdcall, $(args...))
    end
    is_windows() && VERSION >= v"0.6-" && return quote
        ccall(($f,OSQP.osqp), $(esc(:stdcall)), $(args...))
    end
end





function osqp_version()
    return unsafe_string(ccall((:osqp_version, OSQP.osqp), Cstring, ()))
end
