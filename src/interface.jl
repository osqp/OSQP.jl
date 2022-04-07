# Wrapper for the low level functions defined in https://github.com/oxfordcontrol/osqp/blob/master/include/osqp.h

# Ensure compatibility between Julia versions with @gc_preserve
@static if isdefined(Base, :GC)
    import Base.GC: @preserve
else
    macro preserve(args...)
        body = args[end]
        esc(body)
    end
end

"""
    Model()

Initialize OSQP model
"""
mutable struct Model
    solver::Ptr{OSQP.OSQPSolver}
    m::Cc_int
    n::Cc_int
    lcache::Vector{Float64} # to facilitate converting l to use OSQP_INFTY
    ucache::Vector{Float64} # to facilitate converting u to use OSQP_INFTY
	isempty::Bool			# a flag to keep track of the model's setup status
    function Model()
        model = new(C_NULL, 0, 0, Float64[], Float64[], true)
        finalizer(OSQP.clean!, model)
        return model

    end


end

"""
    setup!(model, P, q, A, l, u, settings)

Perform OSQP solver setup of model `model`, using the inputs `P`, `q`, `A`, `l`, `u`.
"""
function setup!(model::OSQP.Model;
        P::Union{SparseMatrixCSC,Nothing} = nothing,
        q::Union{Vector{Float64},Nothing} = nothing,
        A::Union{SparseMatrixCSC,Nothing} = nothing,
        l::Union{Vector{Float64},Nothing} = nothing,
        u::Union{Vector{Float64},Nothing} = nothing,
        settings...)

    # Check problem dimensions
    if P == nothing
        if q != nothing
            n = length(q)
        elseif A != nothing
            n = size(A, 2)
        else
            error("The problem does not have any variables!")
        end

    else
        n = size(P, 1)
    end

    if A == nothing
        m = 0
    else
        m = size(A, 1)
    end


    # Check if parameters are nothing
    if ((A == nothing) & ( (l != nothing) | (u != nothing))) |
        ((A != nothing) & ((l == nothing) & (u == nothing)))
        error("A must be supplied together with l and u")
    end

    if (A != nothing) & (l == nothing)
        l = -Inf * ones(m)
    end
    if (A != nothing) & (u == nothing)
        u = Inf * ones(m)
    end

    if P == nothing
        P = sparse([], [], [], n, n)
    end
    if q == nothing
        q = zeros(n)
    end
    if A == nothing
        A = sparse([], [], [], m, n)
        l = zeros(m)
        u = zeros(m)
    end


    # Check if dimensions are correct
    if length(q) != n
        error("Incorrect dimension of q")
    end
    if length(l) != m
        error("Incorrect dimensions of l")
    end
    if length(u) != m
        error("Incorrect dimensions of u")
    end


    # Check or sparsify matrices
    if !issparse(P)
        @warn("P is not sparse. Sparsifying it now (it might take a while)")
        P = sparse(P)
    end
    if !issparse(A)
        @warn("A is not sparse. Sparsifying it now (it might take a while)")
        A = sparse(A)
    end

    # Constructing upper triangular from P
    if !istriu(P)
        P = triu(P)
    end

    # Convert lower and upper bounds from Julia infinity to OSQP infinity
    u = min.(u, OSQP_INFTY)
    l = max.(l, -OSQP_INFTY)

    # Resize caches
    resize!(model.lcache, m)
    resize!(model.ucache, m)

    # Create managed matrices to avoid segfaults (See SCS.jl)
    managedP = OSQP.ManagedCcsc(P)
    managedA = OSQP.ManagedCcsc(A)

    # Get managed pointers (Ref) Pdata and Adata
    Pdata = Ref(OSQP.Ccsc(managedP))
    Adata = Ref(OSQP.Ccsc(managedA))

    # Create OSQP settings
    settings_dict = Dict{Symbol,Any}()
    if !isempty(settings)
        for (key, value) in settings
            settings_dict[key] = value
        end
    end

    stgs = OSQP.OSQPSettings(settings_dict)

    @preserve managedP Pdata managedA Adata q l u begin

    # # Perform setup
    solver = Ref{Ptr{OSQP.OSQPSolver}}()
    exitflag = ccall(
        (:osqp_setup, OSQP.osqp),
        OSQP.Cc_int,

        (Ptr{Ptr{OSQP.OSQPSolver}}, Ptr{OSQP.Ccsc}, Ptr{Cdouble},
        Ptr{OSQP.Ccsc}, Ptr{Cdouble}, Ptr{Cdouble}, OSQP.Cc_int, OSQP.Cc_int,
        Ptr{OSQP.OSQPSettings}),

        solver,
        Base.unsafe_convert(Ptr{OSQP.Ccsc}, Pdata),
        pointer(q),
        Base.unsafe_convert(Ptr{OSQP.Ccsc}, Adata),
        pointer(l),
        pointer(u),
        m,
        n,
        Ref(stgs)
    )
	model.solver = solver[]

    end

    if exitflag != 0
        error("Error in OSQP setup")
    end

    model.isempty = false
    model.m = m
    model.n = n

end


function solve!(model::OSQP.Model, results::Results = Results())

	model.isempty && throw(ErrorException("You are trying to solve an empty model. Please setup the model before calling solve!()."))
    ccall((:osqp_solve, OSQP.osqp), OSQP.Cc_int,
        (Ref{OSQP.OSQPSolver}, ), model.solver)

    _solver = unsafe_load(model.solver)
    solution = unsafe_load(_solver.solution)
    info = results.info
    copyto!(info, unsafe_load(_solver.info))
    n = model.n
    m = model.m
    resize!(results, n, m)
    has_solution = false
    for status in SOLUTION_PRESENT
        info.status == status && (has_solution = true; break)
    end
    if has_solution
        # If solution exists, copy it
        unsafe_copyto!(pointer(results.x), solution.x, n)
        unsafe_copyto!(pointer(results.y), solution.y, m)
        fill!(results.prim_inf_cert, NaN)
        fill!(results.dual_inf_cert, NaN)
    else
        # else fill with NaN and return certificates of infeasibility
        fill!(results.x, NaN)
        fill!(results.y, NaN)
        if info.status == :Primal_infeasible || info.status == :Primal_infeasible_inaccurate
            unsafe_copyto!(pointer(results.prim_inf_cert), solution.prim_inf_cert, m)
            fill!(results.dual_inf_cert, NaN)
        elseif info.status == :Dual_infeasible || info.status == :Dual_infeasible_inaccurate
            fill!(results.prim_inf_cert, NaN)
            unsafe_copyto!(pointer(results.dual_inf_cert), solution.dual_inf_cert, n)
        else
            fill!(results.prim_inf_cert, NaN)
            fill!(results.dual_inf_cert, NaN)
        end
    end

    if info.status == :Non_convex
        info.obj_val = NaN
    end

    results

end


function version()
    return unsafe_string(ccall((:osqp_version, OSQP.osqp), Cstring, ()))
end

function clean!(model::OSQP.Model)
    exitflag = ccall((:osqp_cleanup, OSQP.osqp), Cc_int,
             (Ptr{OSQP.OSQPSolver},), model.solver)
    if exitflag != 0
        error("Error in OSQP cleanup")
    end
end

function update_q!(model::OSQP.Model, q::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    if length(q) != n
        error("q must have length n = $(n)")
    end
    exitflag = ccall((:osqp_update_data_vec, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, q, C_NULL, C_NULL)
    if exitflag != 0 error("Error updating q") end
end

function update_l!(model::OSQP.Model, l::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    if length(l) != m
        error("l must have length m = $(m)")
    end
    model.lcache .= max.(l, -OSQP_INFTY)
    exitflag = ccall((:osqp_update_data_vec, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, C_NULL, model.lcache, C_NULL)
    if exitflag != 0 error("Error updating l") end
end

function update_u!(model::OSQP.Model, u::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    if length(u) != m
        error("u must have length m = $(m)")
    end
    model.ucache .= min.(u, OSQP_INFTY)
    exitflag = ccall((:osqp_update_data_vec, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, C_NULL, C_NULL, model.ucache)
    if exitflag != 0 error("Error updating u") end
end

function update_bounds!(model::OSQP.Model, l::Vector{Float64}, u::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    if length(l) != m
        error("l must have length m = $(m)")
    end
    if length(u) != m
        error("u must have length m = $(m)")
    end
    model.lcache .= max.(l, -OSQP_INFTY)
    model.ucache .= min.(u, OSQP_INFTY)
    exitflag = ccall((:osqp_update_data_vec, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, C_NULL, model.lcache, model.ucache)
    if exitflag != 0 error("Error updating bounds l and u") end
end

prep_idx_vector_for_ccall(idx::Nothing, n::Int, namesym::Symbol) = C_NULL
function prep_idx_vector_for_ccall(idx::Vector{Int}, n::Int, namesym::Symbol)
    if length(idx) != n
        error("$(namesym) and $(namesym)_idx must have the same length")
    end
    idx .-= 1 # Shift indexing to match C
    idx
end

restore_idx_vector_after_ccall!(idx::Nothing) = nothing
function restore_idx_vector_after_ccall!(idx::Vector{Int})
    idx .+= 1 # Unshift indexing
    nothing
end

function update_P!(model::OSQP.Model, Px::Vector{Float64}, Px_idx::Union{Vector{Int}, Nothing})
    Px_idx_prepped = prep_idx_vector_for_ccall(Px_idx, length(Px), :P)
    exitflag = ccall((:osqp_update_data_mat, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int),
        model.solver, Px, Px_idx_prepped, length(Px), C_NULL, C_NULL, 0)
    restore_idx_vector_after_ccall!(Px_idx)
    if exitflag != 0 error("Error updating P") end
end

function update_A!(model::OSQP.Model, Ax::Vector{Float64}, Ax_idx::Union{Vector{Int}, Nothing})
    Ax_idx_prepped = prep_idx_vector_for_ccall(Ax_idx, length(Ax), :A)
    exitflag = ccall((:osqp_update_data_mat, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int),
        model.solver, C_NULL, C_NULL, 0, Ax, Ax_idx_prepped, length(Ax))
    restore_idx_vector_after_ccall!(Ax_idx)
    if exitflag != 0 error("Error updating A") end
end

function update_P_A!(model::OSQP.Model, Px::Vector{Float64}, Px_idx::Union{Vector{Int}, Nothing}, Ax::Vector{Float64}, Ax_idx::Union{Vector{Int}, Nothing})
    Px_idx_prepped = prep_idx_vector_for_ccall(Px_idx, length(Px), :P)
    Ax_idx_prepped = prep_idx_vector_for_ccall(Ax_idx, length(Ax), :A)
    exitflag = ccall((:osqp_update_data_mat, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int, Ptr{Cdouble}, Ptr{Cc_int}, Cc_int),
        model.solver, Px, Px_idx_prepped, length(Px), Ax, Ax_idx_prepped, length(Ax))
    restore_idx_vector_after_ccall!(Ax_idx)
    restore_idx_vector_after_ccall!(Px_idx)
    if exitflag != 0 error("Error updating P and A") end
end

function update!(model::OSQP.Model; q = nothing, l = nothing, u = nothing, Px = nothing, Px_idx = nothing, Ax = nothing, Ax_idx = nothing)
    # q
    if q != nothing
        update_q!(model, q)
    end

    # l and u
    if l != nothing && u != nothing
        update_bounds!(model, l, u)
    elseif l != nothing
        update_l!(model, l)
    elseif u != nothing
        update_u!(model, u)
    end

    # P and A
    if Px != nothing && Ax != nothing
        update_P_A!(model, Px, Px_idx, Ax, Ax_idx)
    elseif Px != nothing
        update_P!(model, Px, Px_idx)
    elseif Ax != nothing
        update_A!(model, Ax, Ax_idx)
    end
end



function update_settings!(model::OSQP.Model; kwargs...)

    new_settings = Dict(kwargs)
    if :rho in keys(new_settings)
        rho = pop!(new_settings, :rho)
        exitflag = ccall((:osqp_update_rho, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Cdouble), model.solver, rho)
        if exitflag != 0 error("Error updating rho") end
    end

    if !isempty(new_settings)
        for (key, _) in new_settings
            if !(key in UPDATABLE_SETTINGS)
                error("$(key) cannot be updated or is not recognized. Not updating anything.")
            end
        end

        data = Dict{Symbol,Any}()
        old_settings = unsafe_load(unsafe_load(model.solver).data)
        for fieldname in fieldnames(OSQP.OSQPSettings)
            data[fieldname] = getfield(old_settings, fieldname)
        end
        data = merge(data, new_settings)

        settings = OSQPSettings(data)
        exitflag = ccall((:osqp_update_settings, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{OSQP.OSQPSettings}), model.solver, Ref(settings))
        if exitflag != 0 error("Error updating settings") end
    end
end

function warm_start_x!(model::OSQP.Model, x::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    length(x) == n || error("Wrong dimension for variable x")
    exitflag = ccall((:osqp_warm_start, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, x, C_NULL)
    exitflag == 0  || error("Error in warm starting x")
    nothing
end

function warm_start_y!(model::OSQP.Model, y::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    length(y) == m || error("Wrong dimension for variable y")
    exitflag = ccall((:osqp_warm_start, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, C_NULL, y)
    exitflag == 0 || error("Error in warm starting y")
    nothing
end

function warm_start_x_y!(model::OSQP.Model, x::Vector{Float64}, y::Vector{Float64})
    (n, m) = OSQP.dimensions(model)
    length(x) == n || error("Wrong dimension for variable x")
    length(y) == m || error("Wrong dimension for variable y")
    exitflag = ccall((:osqp_warm_start, OSQP.osqp), Cc_int, (Ptr{OSQP.OSQPSolver}, Ptr{Cdouble}, Ptr{Cdouble}), model.solver, x, y)
    exitflag == 0 || error("Error in warm starting x and y")
    nothing
end


function warm_start!(model::OSQP.Model; x::Union{Vector{Float64}, Nothing} = nothing, y::Union{Vector{Float64}, Nothing} = nothing)
    if x isa Vector{Float64} && y isa Vector{Float64}
        warm_start_x_y!(model, x, y)
    elseif x isa Vector{Float64}
        warm_start_x!(model, x)
    elseif y isa Vector{Float64}
        warm_start_y!(model, y)
    end
end

# Auxiliary low-level functions
"""
    dimensions(model::OSQP.Model)

Obtain problem dimensions from OSQP model
"""
function dimensions(model::OSQP.Model)
    # TODO: Why was this function returning (n, m) instead of (m, n)?
    return model.n, model.m
end


function linsys_solver_str_to_int!(settings_dict::Dict{Symbol,Any})
    linsys_str = get(settings_dict, :linsys_solver, nothing)

    if linsys_str != nothing
         # Check type
        if !isa(linsys_str, String)
            # If linsys_str is already an expected value, do nothing
            if !(linsys_str in (OSQP_DIRECT_SOLVER, OSQP_INDIRECT_SOLVER))
                error("Linear system solver not recognized.")
            else
                return
            end
        end

         # Convert to lower case
        linsys_str = lowercase(linsys_str)

        if linsys_str == "direct"
            settings_dict[:linsys_solver] = OSQP_DIRECT_SOLVER
        elseif linsys_str == "indirect"
            settings_dict[:linsys_solver] = OSQP_INDIRECT_SOLVER
        else
            error("Linear system solver not recognized.")
        end
    end
end