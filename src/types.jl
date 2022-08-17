# Types defined in types.h
# https://github.com/oxfordcontrol/osqp/blob/master/include/types.h

const CFloats = Union{Cfloat, Cdouble}
const CInts = Union{Cint, Clonglong}

struct Ccsc{FT<:CFloats, IT<:CInts}
    m::IT
    n::IT
    p::Ptr{IT}
    i::Ptr{IT}
    x::Ptr{FT}
    nzmax::IT
    nz::IT
end

struct ManagedCcsc{FT<:CFloats, IT<:CInts}
    m::IT
    n::IT
    p::Vector{IT}
    i::Vector{IT}
    x::Vector{FT}
    nzmax::IT
    nz::IT
end

# Construct ManagedCcsc matrix from SparseMatrixCSC
function ManagedCcsc(M::SparseMatrixCSC, FT, IT)

    # Get dimensions
    m = M.m
    n = M.n

    # Get vectors of data, rows indices and column pointers
    x = convert(Array{FT,1}, M.nzval)
    # C is 0 indexed
    i = convert(Array{IT,1}, M.rowval .- 1)
    # C is 0 indexed
    p = convert(Array{IT,1}, M.colptr .- 1)

    # Create new ManagedCcsc matrix
    return ManagedCcsc{FT,IT}(m, n, p, i, x, length(M.nzval), -1)
end

function Base.convert(::Type{SparseMatrixCSC}, c::OSQP.Ccsc)
    m = c.m
    n = c.n
    nzmax = c.nzmax
    nzval = [unsafe_load(c.x, i) for i in 1:nzmax]
    rowval = [unsafe_load(c.i, i) for i in 1:nzmax] .+ 1
    colptr = [unsafe_load(c.p, i) for i in 1:(n+1)] .+ 1
    return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

# Returns an Ccsc matrix. The vectors are *not* GC tracked in the struct.
# Use this only when you know that the managed matrix will outlive the Ccsc
# matrix.
function Ccsc(m::ManagedCcsc{FT,IT}) where {FT<:CFloats, IT<:CInts}
    return Ccsc{FT,IT}(
        m.m,
        m.n,
        pointer(m.p),
        pointer(m.i),
        pointer(m.x),
        m.nzmax,
        m.nz
    )
end

struct OSQPSolution{FT<:CFloats, IT<:CInts}
    x::Ptr{FT}
    y::Ptr{FT}
    prim_inf_cert::Ptr{FT}
    dual_inf_cert::Ptr{FT}
end

mutable struct OSQPSettings{FT<:CFloats, IT<:CInts}
    device::IT
    linsys_solver::Cint  # Enum type
    verbose::IT
    warm_starting::IT
    scaling::IT
    polishing::IT
    rho::FT
    rho_is_vec::IT
    sigma::FT
    alpha::FT
    cg_max_iter::IT
    cg_tol_reduction::IT
    cg_tol_fraction::FT
    adaptive_rho::IT
    adaptive_rho_interval::IT
    adaptive_rho_fraction::FT
    adaptive_rho_tolerance::FT
    max_iter::IT
    eps_abs::FT
    eps_rel::FT
    eps_prim_inf::FT
    eps_dual_inf::FT
    scaled_termination::IT
    check_termination::IT
    time_limit::FT
    delta::FT
    polish_refine_iter::IT

    function OSQPSettings{FT,IT}() where {FT, IT}
        return new{FT, IT}()
    end
end

# UPDATABLE_SETTINGS
# Since OSQP's osqp_update_settings function takes a pointer to an OSQPSettings object,
# all fields of the OSQPSettings struct are updatable.
const UPDATABLE_SETTINGS = fieldnames(OSQPSettings)

function OSQPSettings(algebra::alg) where {alg<:OSQPAlgebra}
    return get_default_settings(algebra)
end

function OSQPSettings(settings_dict::Dict{Symbol,Any}, algebra::alg = OSQPBuiltinAlgebra()) where {alg<:OSQPAlgebra}
    settings = get_default_settings(algebra)

    # Convert linsys_solver string to number
    linsys_solver_str_to_int!(settings_dict)

    for (key, val) in settings_dict
        setfield!(settings, key, convert(fieldtype(typeof(settings), key), val))
    end

    return settings
end

struct OSQPInfo{FT<:CFloats, IT<:CInts}
    # We need to allocate 32 bytes for a character string, so we allocate 256 bits
    # of integer instead
    # TODO: Find a better way to do this
    status::NTuple{32,Cchar}
    status_val::IT
    status_polish::IT
    obj_val::FT
    prim_res::FT
    dual_res::FT
    iter::IT
    rho_updates::IT
    rho_estimate::FT
    setup_time::FT
    solve_time::FT
    update_time::FT
    polish_time::FT
    run_time::FT
end

# OSQPWorkspace is an internal struct and remains opaque to us
mutable struct OSQPWorkspace
end

mutable struct OSQPSolver{FT<:CFloats, IT<:CInts}
    data::Ptr{OSQP.OSQPSettings{FT,IT}}
    solution::Ptr{OSQP.OSQPSolution{FT,IT}}
    info::Ptr{OSQP.OSQPInfo{FT,IT}}
    work::Ptr{OSQP.OSQPWorkspace}
end

mutable struct Info{FT,IT}
    status::Symbol
    status_val::IT
    status_polish::IT
    obj_val::FT
    prim_res::FT
    dual_res::FT
    iter::IT
    rho_updates::IT
    rho_estimate::FT
    setup_time::FT
    solve_time::FT
    update_time::FT
    polish_time::FT
    run_time::FT

    function Info{FT,IT}() where {FT,IT}
        return new()
    end
end

function copyto!(info::Info, cinfo::OSQPInfo{FT,IT}) where {FT<:CFloats, IT<:CInts}
    info.status = OSQP.status_map[cinfo.status_val]
    info.status_val = cinfo.status_val
    info.status_polish = cinfo.status_polish
    info.obj_val = cinfo.obj_val
    info.prim_res = cinfo.prim_res
    info.dual_res = cinfo.dual_res
    info.iter = cinfo.iter
    info.rho_updates = cinfo.rho_updates
    info.rho_estimate = cinfo.rho_estimate
    info.setup_time = cinfo.setup_time
    info.solve_time = cinfo.solve_time
    info.update_time = cinfo.update_time
    info.polish_time = cinfo.polish_time
    info.run_time = cinfo.run_time
    return info
end

mutable struct Results{FT,IT}
    x::Vector{FT}
    y::Vector{FT}
    info::OSQP.Info{FT,IT}
    prim_inf_cert::Vector{FT}
    dual_inf_cert::Vector{FT}

    function Results{FT,IT}() where {FT,IT}
        res = new{FT,IT}()
        res.info = Info{FT,IT}()
        res.x = Vector{FT}(undef,0)
        res.y = Vector{FT}(undef,0)
        res.prim_inf_cert = Vector{FT}(undef,0)
        res.dual_inf_cert = Vector{FT}(undef,0)

        return res
    end
end

function Base.resize!(results::Results, n::IN, m::IN) where {IN <: Union{Int32, Int64}}
    resize!(results.x, n)
    resize!(results.y, m)
    resize!(results.prim_inf_cert, m)
    resize!(results.dual_inf_cert, n)
    return results
end
