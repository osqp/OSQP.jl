using Base.Meta: isexpr, quot
using MacroTools: postwalk

"""
    OSQPAlgebra{FT,IT}

The supertype for all backend algebra library implementations of OSQP to implement
a library pointer type underneath. The type parameters `FT` and `IT` refer to the
types for the floating-point numbers and integers the underlying library uses,
respectively.
"""
abstract type OSQPAlgebra{FT,IT} end


"""
    _osqp_jlfunc_name(func)

Create the name that will be given to the Julia function to wrap the ccall of `func`.
"""
_osqp_jlfunc_name(func) = Symbol("_ccall_" * replace(string(func), r":" => ""))


"""
    cfunc

Internal type used to represent a C function and its return type/arguments.
"""
struct cfunc
    # The function name
    name::Symbol

    # The function return type
    rettype::Union{Symbol,Expr}

    # The names of the arguments to the function
    args::Vector{Symbol}

    # The types of the arguments to the function
    argt::Vector{Union{Symbol,Expr}}
end

# A dictionary that will contain information on all C functions exposed in the OSQP library
# This is automatically added to by the @cprototype macro.
const osqp_api = Dict{String,cfunc}()


"""
    @cprototype sig

Define a prototype for the C function given by `sig`.

When defining the function, `Tint` and `Tfloat` should be used to represent the types
`c_float` and `c_int` in the underlying C API, respectively. These will be automatically
converted to the proper types when the definition of the function is created using [`@cdefinition`](@ref).
"""
macro cprototype(sig)
    expr = _cprototype(__module__, __source__, sig)
    return esc(expr)
end

function _cprototype(__module__, __source__, sig::Expr)
    sig.head === :(::) || error("return type required on function signature")

    # Pull apart return-type and rest of function declaration
    rettype = sig.args[2]::Union{Symbol,Expr}
    funcsig = sig.args[1]
    isexpr(funcsig, :call) || error("expected function-like expression, found `", funcsig, "`")
    funcsig = funcsig::Expr

    # Extract function name and argument list
    jlfuncname = funcsig.args[1]::Symbol
    funcargs = funcsig.args[2:end]

    # Pull apart argument names and types
    args = Vector{Symbol}()
    argt = Vector{Union{Expr,Symbol}}()
    for ii in 1:length(funcargs)
        argex = funcargs[ii]
        if !isexpr(argex, :(::)) || !(argex.args[1] isa Symbol)
            error("expected `name::type` expression in argument ", ii, ", got ", funcargs[ii])
        end
        push!(args, argex.args[1])
        push!(argt, argex.args[2])
    end

    # Store the function in the function database
    sfname = String(jlfuncname)
    osqp_api[sfname] = cfunc(jlfuncname, rettype, args, argt)

    # Create a base method that provides an error and something for the actual algebra backends
    # to extend with their own type
    fname = _osqp_jlfunc_name(jlfuncname)

    f = quote
            function $(fname)(algebra::alg, $(args...)) where {alg <: OSQPAlgebra}
                error(":" * $(sfname) * " is not implemented for the " * string(typeof(algebra)) * " algebra library")
            end
        end

    return f
end

"""
    @cdefinition backend, library, Tfloat, Tint, funcname

Define a method that will call the C function `funcname` (which was previously defined
by [`@cprototype`](@ref) in the main OSQP package) inside the library with the handle
given by `library` and associated it with the linear algebra backend `backend`.

The Tsub dictionary specified a mapping from symbol => symbol for the types specified in the
prototype to be replaced with a given concrete type.
"""
macro cdefinition(backend, library, funcname, Tsub)
    expr = _cdefinition(__module__, __source__, backend, library, funcname, Tsub)
    return esc(expr)
end

function _cdefinition(__module__, __source__, backend, library, funcname, Tsub)
    sfname = String(funcname)

    Tdict = __module__.eval(Tsub)

    # Modify the stored API definitions to have the proper types for the backend library
    # by replacing Tint and Tfloat.
    rawfunc = OSQP.osqp_api[sfname]

    # The types can get pretty complicated, so we need to ensure we recurse into all expressions
    # to get all the uses of Tfloat and Tint, hence the use of the postwalk function.
    rep_argt    = rawfunc.argt
    rep_rettype = rawfunc.rettype

    for (Tholder, Tactual) in Tdict
        rep_argt    = map(x -> postwalk(x -> x == Tholder ? Tactual : x, x), rep_argt)
        rep_rettype = postwalk(x -> x == Tholder ? Tactual : x, rep_rettype)
    end

    repfunc = cfunc(rawfunc.name, rep_rettype, rawfunc.args, rep_argt)

    # Build the ccall itself
    cfunclib = Expr(:tuple, quot(repfunc.name), library)
    ccallexpr = :(ccall($(cfunclib), $(repfunc.rettype), ($(repfunc.argt...),), $(repfunc.args...)))

    # Assemble the full function
    jlfuncsig = Expr(:call, _osqp_jlfunc_name(repfunc.name), :(::$(backend)), repfunc.args...)
    jlfuncbody = Expr(:block, __source__, :(return $ccallexpr))
    jlfuncexpr = Expr(:function, jlfuncsig, jlfuncbody)

    # The function extends an OSQP function, so there must be an import statement before defining the function,
    # but only if we aren't actually in the OSQP module already.
    if __module__ != OSQP
        jlfuncimport = :(import OSQP: $(_osqp_jlfunc_name(repfunc.name)))
        jlfullexpr = Expr(:block, __source__, jlfuncimport, jlfuncexpr)
        return jlfullexpr
    end

    return jlfuncexpr
end

"""
    @osqp_ccall(func, alg::algebra, args...) where {backend <: OSQPAlgebra}

Call a function `func` inside the OSQP backend linear algebra library specified by `alg` and
pass `args...` to the underlying C function.

All functions called using this macro must have previously been defined using the [`@cprototype`](@ref)
in the main OSQP package and the [`@cdefinition`](@ref) in any backend packages providing it.
"""
macro osqp_ccall(func, alg, args...)
    fname = _osqp_jlfunc_name(func)

    return esc(:($(fname)($(alg), $(args...))))
end
