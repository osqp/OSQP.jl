# Types defined in types.h
# https://github.com/oxfordcontrol/osqp/blob/master/include/types.h

# Integer type from C
if Sys.WORD_SIZE == 64   # 64bit system
    Cc_int = Clonglong
else  # 32bit system
    Cc_int = Cint
end

struct Ccsc
    nzmax::Cc_int
    m::Cc_int
    n::Cc_int
    p::Ptr{Cc_int}
    i::Ptr{Cc_int}
    x::Ptr{Cdouble}
    nz::Cc_int
end


struct ManagedCcsc
    nzmax::Cc_int
    m::Cc_int
    n::Cc_int
    p::Vector{Cc_int}
    i::Vector{Cc_int}
    x::Vector{Cdouble}
    nz::Cc_int

end

# Construct ManagedCcsc matrix from SparseMatrixCSC
function ManagedCcsc(M::SparseMatrixCSC)

    # Get dimensions
    m = M.m
    n = M.n

    # Get vectors of data, rows indices and column pointers
    x = convert(Array{Float64, 1}, M.nzval)
    # C is 0 indexed
    i = convert(Array{Cc_int, 1}, M.rowval .- 1)
    # C is 0 indexed
    p = convert(Array{Cc_int, 1}, M.colptr .- 1)
    
    # Create new ManagedCcsc matrix
    ManagedCcsc(length(M.nzval), m, n, p, i, x, -1)
end

# Returns an Ccsc matrix. The vectors are *not* GC tracked in the struct.
# Use this only when you know that the managed matrix will outlive the Ccsc
# matrix.
Ccsc(m::ManagedCcsc) =
    Ccsc(m.nzmax, m.m, m.n, pointer(m.p), pointer(m.i), pointer(m.x), m.nz)


struct Solution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
end

# Internal C type for info
# N.B. This is not the one returned to the user!
struct CInfo
    iter::Cc_int
    # We need to allocate 32 bytes for a character string, so we allocate 256 bits
    # of integer instead
    # TODO: Find a better way to do this
    status::NTuple{32, Cchar}
    status_val::Cc_int
    status_polish::Cc_int
    obj_val::Cdouble
    pri_res::Cdouble
    dua_res::Cdouble
    setup_time::Cdouble
    solve_time::Cdouble
    polish_time::Cdouble
    run_time::Cdouble
    rho_updates::Cc_int
    rho_estimate::Cdouble
end

struct Data
    n::Cc_int
    m::Cc_int
    P::Ptr{Ccsc}
    A::Ptr{Ccsc}
    q::Ptr{Cdouble}
    l::Ptr{Cdouble}
    u::Ptr{Cdouble}
end


struct Settings
    rho::Cdouble
    sigma::Cdouble
    scaling::Cc_int
    adaptive_rho::Cc_int
    adaptive_rho_interval::Cc_int
    adaptive_rho_tolerance::Cdouble
    adaptive_rho_fraction::Cdouble
    max_iter::Cc_int
    eps_abs::Cdouble
    eps_rel::Cdouble
    eps_prim_inf::Cdouble
    eps_dual_inf::Cdouble
    alpha::Cdouble
    linsys_solver::Cint  # Enum type
    delta::Cdouble
    polish::Cc_int
    polish_refine_iter::Cc_int
    verbose::Cc_int
    scaled_termination::Cc_int
    check_termination::Cc_int
    warm_start::Cc_int
end

function Settings()
    s = Ref{OSQP.Settings}()
    ccall((:set_default_settings, OSQP.osqp), Nothing,
          (Ref{OSQP.Settings},), s)
    return s[]
end

function Settings(settings_dict::Dict{Symbol, Any})
#  function Settings(settings::Base.Iterators.IndexValue)
#  function Settings(settings::Array{Any, 1})
    default_settings = OSQP.Settings()


       # Convert linsys_solver string to number
       linsys_solver_str_to_int!(settings_dict)

    # Get list with elements of default and user settings
    # If setting is in the passed settings (settings_dict), 
    # then convert type to the right type. Otherwise just take
    # the default one
    settings_list = [setting in keys(settings_dict) ?
             convert(fieldtype(typeof(default_settings), setting), settings_dict[setting]) :
             getfield(default_settings, setting)
             for setting in fieldnames(typeof(default_settings))]

    # Create new settings with new dictionary
    s = OSQP.Settings(settings_list...)
    return s

end



struct Workspace
    data::Ptr{OSQP.Data}
    linsys_solver::Ptr{Nothing}
    pol::Ptr{Nothing}

    rho_vec::Ptr{Cdouble}
    rho_inv_vec::Ptr{Cdouble}
    constr_type::Ptr{Cc_int}

    # Iterates
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    z::Ptr{Cdouble}
    xz_tilde::Ptr{Cdouble}
    x_prev::Ptr{Cdouble}
    z_prev::Ptr{Cdouble}

    # Primal and dual residuals
    Ax::Ptr{Cdouble}
    Px::Ptr{Cdouble}
    Aty::Ptr{Cdouble}

    # Primal infeasibility
    delta_y::Ptr{Cdouble}
    Atdelta_y::Ptr{Cdouble}

    # Dual infeasibility
    delta_x::Ptr{Cdouble}
    Pdelta_x::Ptr{Cdouble}
    Adelta_x::Ptr{Cdouble}


    # Scaling
    D_temp::Ptr{Cdouble}
    D_temp_A::Ptr{Cdouble}
    E_temp::Ptr{Cdouble}

    settings::Ptr{OSQP.Settings}
    scaling::Ptr{Nothing}
    solution::Ptr{OSQP.Solution}
    info::Ptr{OSQP.CInfo}

    timer::Ptr{Nothing}
    first_run::Cc_int
    summary_printed::Cc_int

end


struct Info
    iter::Int64
    status::Symbol
    status_val::Int64
    status_polish::Int64
    obj_val::Float64
    pri_res::Float64
    dua_res::Float64
    setup_time::Float64
    solve_time::Float64
    polish_time::Float64
    run_time::Float64
    rho_updates::Int64
    rho_estimate::Float64

    function Info(cinfo::CInfo)
        status = OSQP.status_map[cinfo.status_val]
        return new(cinfo.iter,
                   status,
                   cinfo.status_val,
                   cinfo.status_polish,
                   cinfo.obj_val,
                   cinfo.pri_res,
                   cinfo.dua_res,
                   cinfo.setup_time,
                   cinfo.solve_time,
                   cinfo.polish_time,
                   cinfo.run_time,
               cinfo.rho_updates,
               cinfo.rho_estimate)
    end
end

struct Results
    x::Vector{Float64}
    y::Vector{Float64}
    info::OSQP.Info
    prim_inf_cert::Union{Vector{Float64}, Nothing}
    dual_inf_cert::Union{Vector{Float64}, Nothing}

end
Results(x, y, info) = Results(x, y, info, nothing, nothing) 



