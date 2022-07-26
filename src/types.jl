# Types defined in types.h
# https://github.com/oxfordcontrol/osqp/blob/master/include/types.h

# Integer type from C
if Sys.WORD_SIZE == 64   # 64bit system
    const Cc_int = Clonglong
else  # 32bit system
    const Cc_int = Cint
end

struct Ccsc
    m::Cc_int
    n::Cc_int
    p::Ptr{Cc_int}
    i::Ptr{Cc_int}
    x::Ptr{Cdouble}
    nzmax::Cc_int
    nz::Cc_int
end

struct ManagedCcsc
    m::Cc_int
    n::Cc_int
    p::Vector{Cc_int}
    i::Vector{Cc_int}
    x::Vector{Cdouble}
    nzmax::Cc_int
    nz::Cc_int
end

# Construct ManagedCcsc matrix from SparseMatrixCSC
function ManagedCcsc(M::SparseMatrixCSC)

    # Get dimensions
    m = M.m
    n = M.n

    # Get vectors of data, rows indices and column pointers
    x = convert(Array{Float64,1}, M.nzval)
    # C is 0 indexed
    i = convert(Array{Cc_int,1}, M.rowval .- 1)
    # C is 0 indexed
    p = convert(Array{Cc_int,1}, M.colptr .- 1)

    # Create new ManagedCcsc matrix
    return ManagedCcsc(m, n, p, i, x, length(M.nzval), -1)
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
function Ccsc(m::ManagedCcsc)
    return Ccsc(
        m.m,
        m.n,
        pointer(m.p),
        pointer(m.i),
        pointer(m.x),
        m.nzmax,
        m.nz
    )
end

struct OSQPSolution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    prim_inf_cert::Ptr{Cdouble}
    dual_inf_cert::Ptr{Cdouble}
end

struct OSQPSettings
    device::Cc_int
    linsys_solver::Cint  # Enum type
    verbose::Cc_int
    warm_starting::Cc_int
    scaling::Cc_int
    polishing::Cc_int
    rho::Cdouble
    rho_is_vec::Cc_int
    sigma::Cdouble
    alpha::Cdouble
    cg_max_iter::Cc_int
    cg_tol_reduction::Cc_int
    cg_tol_fraction::Cdouble
    adaptive_rho::Cc_int
    adaptive_rho_interval::Cc_int
    adaptive_rho_fraction::Cdouble
    adaptive_rho_tolerance::Cdouble
    max_iter::Cc_int
    eps_abs::Cdouble
    eps_rel::Cdouble
    eps_prim_inf::Cdouble
    eps_dual_inf::Cdouble
    scaled_termination::Cc_int
    check_termination::Cc_int
    time_limit::Cdouble
    delta::Cdouble
    polish_refine_iter::Cc_int
end

# UPDATABLE_SETTINGS
# Since OSQP's osqp_update_settings function takes a pointer to an OSQPSettings object,
# all fields of the OSQPSettings struct are updatable.
const UPDATABLE_SETTINGS = fieldnames(OSQPSettings)

function OSQPSettings(algebra::alg) where {alg<:OSQPAlgebra}
    return get_default_settings(algebra)
end

function OSQPSettings(settings_dict::Dict{Symbol,Any}, algebra::alg = OSQPBuiltinAlgebra()) where {alg<:OSQPAlgebra}
    default_settings = get_default_settings(algebra)

    # Convert linsys_solver string to number
    linsys_solver_str_to_int!(settings_dict)

    # Get list with elements of default and user settings
    # If setting is in the passed settings (settings_dict),
    # then convert type to the right type. Otherwise just take
    # the default one
    settings_list = [
        setting in keys(settings_dict) ?
        convert(
            fieldtype(typeof(default_settings), setting),
            settings_dict[setting],
        ) : getfield(default_settings, setting) for
        setting in fieldnames(typeof(default_settings))
    ]

    # Create new settings with new dictionary
    s = OSQP.OSQPSettings(settings_list...)
    return s
end

struct OSQPInfo
    # We need to allocate 32 bytes for a character string, so we allocate 256 bits
    # of integer instead
    # TODO: Find a better way to do this
    status::NTuple{32,Cchar}
    status_val::Cc_int
    status_polish::Cc_int
    obj_val::Cdouble
    prim_res::Cdouble
    dual_res::Cdouble
    iter::Cc_int
    rho_updates::Cc_int
    rho_estimate::Cdouble
    setup_time::Cdouble
    solve_time::Cdouble
    update_time::Cdouble
    polish_time::Cdouble
    run_time::Cdouble
end

# OSQPWorkspace is an internal struct and remains opaque to us
mutable struct OSQPWorkspace
end

mutable struct OSQPSolver
    data::Ptr{OSQP.OSQPSettings}
    solution::Ptr{OSQP.OSQPSolution}
    info::Ptr{OSQP.OSQPInfo}
    work::Ptr{OSQP.OSQPWorkspace}
end

mutable struct Info
    status::Symbol
    status_val::Int64
    status_polish::Int64
    obj_val::Float64
    prim_res::Float64
    dual_res::Float64
    iter::Int64
    rho_updates::Int64
    rho_estimate::Float64
    setup_time::Float64
    solve_time::Float64
    update_time::Float64
    polish_time::Float64
    run_time::Float64

    Info() = new()
end

function copyto!(info::Info, cinfo::OSQPInfo)
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

mutable struct Results
    x::Vector{Float64}
    y::Vector{Float64}
    info::OSQP.Info
    prim_inf_cert::Vector{Float64}
    dual_inf_cert::Vector{Float64}
end

Results() = Results(Float64[], Float64[], Info(), Float64[], Float64[])

function Base.resize!(results::Results, n::Int, m::Int)
    resize!(results.x, n)
    resize!(results.y, m)
    resize!(results.prim_inf_cert, m)
    resize!(results.dual_inf_cert, n)
    return results
end
