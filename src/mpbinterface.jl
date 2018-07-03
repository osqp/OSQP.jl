module OSQPMathProgBaseInterface

using Compat
using Compat.LinearAlgebra
using Compat.SparseArrays
using OSQP: Model, Results, setup!, solve!, update!, clean!, update_settings!, warm_start!
using MathProgBase

@static if VERSION < v"0.7-"
    blockdiag(A...) = blkdiag(A...)
end

struct OSQPSolver <: MathProgBase.AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

OSQPSolver(; kwargs...) = OSQPSolver(Dict{Symbol, Any}(k => v for (k, v) in kwargs))

mutable struct OSQPMathProgModel <: MathProgBase.AbstractLinearQuadraticModel
    # Problem settings
    settings::Dict{Symbol, Any}

    # Inner OSQP model structure
    inner::Model

    # Setup needed
    perform_setup::Bool

    # Problem data in the form
    #
    # minimize    (1/2) x' P x + q' x
    # subject to  lb <= A x <= ub
    #             l <= x <= u
    sense::Symbol
    P::SparseMatrixCSC{Float64, Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    lb::Vector{Float64}
    ub::Vector{Float64}
    l::Vector{Float64}
    u::Vector{Float64}

    # Warm start auxiliary variables
    x_ws::Vector{Float64}  # Warm start vector
    # Call warm start function in optimize!.
    # NB. Warmstart is enabled internally by default when changing problem data
    run_warmstart::Bool

    # Solved
    solved::Bool

    # Results
    results::Results

    function OSQPMathProgModel(settings::AbstractDict{Symbol, Any})
        P = spzeros(Float64, 0, 0)
        q = Float64[]
        A = spzeros(Float64, 0, 0)
        lb = Float64[]
        ub = Float64[]
        l = Float64[]
        u = Float64[]
        perform_setup = true
        x_ws = Float64[]
        run_warmstart = false
        solved = false
        new(settings, Model(), perform_setup,
            :Min, P, q, A, lb, ub, l, u,
            x_ws, run_warmstart, solved) # leave other fields undefined
    end
end

function MathProgBase.numvar(model::OSQPMathProgModel)
    # Check problem dimensions
    if size(model.P, 2) == 0
        if size(model.A, 2) != 0
            n = size(model.A, 2)
        else
            error("the model has no variables")
        end
    else
        n = size(model.P, 1)
    end

    return n
end

MathProgBase.numconstr(model::OSQPMathProgModel) = size(model.A, 1)

# Initialize
function initialize_dimensions!(model::OSQPMathProgModel, n, m)
    # Initialize empty vectors and matrices with the right dimensions
    model.P = spzeros(n, n)
    model.q = zeros(n)
    model.A = spzeros(m, n)
    model.lb = zeros(m)
    model.ub = zeros(m)
    model.l = zeros(n)
    model.u = zeros(n)
    model.x_ws = zeros(n)

    model
end


# Maybe add resetproblem function to remove results and reset perform_setup?
checksolved(model::OSQPMathProgModel) = model.solved || error("Model has not been solved.")

# Reset problem when setup has to be performed again
resetproblem(model::OSQPMathProgModel) = (model.solved = false)  # Clear model.results!

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
MathProgBase.LinearQuadraticModel(solver::OSQPSolver) = OSQPMathProgModel(solver.settings)
MathProgBase.getsolution(model::OSQPMathProgModel) = (checksolved(model); model.results.x)
function MathProgBase.getobjval(model::OSQPMathProgModel)
    checksolved(model)

    if model.sense == :Min
        return model.results.info.obj_val
    else
        return -model.results.info.obj_val
    end
end

function MathProgBase.optimize!(model::OSQPMathProgModel)
    # Perform setup only if necessary
    if model.perform_setup
        n = MathProgBase.numvar(model)
        setup!(model.inner; P = model.P, q = model.q,
               A = [model.A; sparse(I, n, n)],
               l = [model.lb; model.l],
               u = [model.ub; model.u],
               model.settings...)
        model.perform_setup = false  # No longer needed to perform setup
    end

    # Call internal warm start function if needed
    if model.run_warmstart
        warm_start!(model.inner, x=model.x_ws)
    end

    # Solve the problem and store the results
    model.results = solve!(model.inner)

    # Problem is now solved
    model.solved = true
end

function MathProgBase.status(model::OSQPMathProgModel)::Symbol
    checksolved(model)
    osqpstatus = model.results.info.status
    status = osqpstatus # if OSQP status can't be mapped to a standard return value, just return as is

    # map to standard return status values:
    osqpstatus == :Solved && (status = :Optimal)
    osqpstatus == :Max_iter_reached && (status = :UserLimit)
    osqpstatus == :Interrupted && (status = :UserLimit) # following Gurobi.jl
    osqpstatus == :Primal_infeasible && (status = :Infeasible)
    osqpstatus == :Primal_infeasible_inaccurate && (status = :Infeasible)
    osqpstatus == :Dual_infeasible && (status = :Unbounded)
    osqpstatus == :Dual_infeasible_inaccurate && (status = :Unbounded)
    osqpstatus == :Time_limit_reached && (status = :UserLimit)

    return status
end

# TODO: getobjbound
# TODO: getobjgap

MathProgBase.getrawsolver(model::OSQPMathProgModel) = model.inner
MathProgBase.getsolvetime(model::OSQPMathProgModel) = (checksolved(model); model.results.info.run_time)


function MathProgBase.setsense!(model::OSQPMathProgModel, sense::Symbol)
    if sense in [:Min, :Max]
        if (sense == :Max) & (model.sense != :Max)
            model.q .= .-model.q
            model.P .= .-model.P

            # Update problem data
            if !model.perform_setup
                update!(model.inner, q=model.q)
                update!(model.inner, Px=triu(model.P).nzval)
            end
        end
        model.sense = sense
    else
        error("sense not recognized")
    end

    resetproblem(model)
end

MathProgBase.getsense(model::OSQPMathProgModel) = model.sense
MathProgBase.freemodel!(model::OSQPMathProgModel) = clean!(model.inner)

function Base.copy(model::OSQPMathProgModel)
    ret = OSQPMathProgModel(model.settings)
    initialize_dimensions!(ret, MathProgBase.numvar(model), MathProgBase.numconstr(model))
    ret.sense = model.sense
    copyto!(ret.P, model.P)
    copyto!(ret.q, model.q)
    copyto!(ret.A, model.A)
    copyto!(ret.lb, model.lb)
    copyto!(ret.ub, model.ub)
    copyto!(ret.l, model.l)
    copyto!(ret.u, model.u)
    ret
end

MathProgBase.setvartype!(model::OSQPMathProgModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
MathProgBase.getvartype(model::OSQPMathProgModel) = fill(:Cont, MathProgBase.numvar(model))

# Set solver independent parameters: Silent is the only supported one for now
function MathProgBase.setparameters!(x::OSQPMathProgModel; Silent = nothing, TimeLimit = nothing)
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent

        # Update silent setting if setup has already been performed
        if !x.perform_setup
            update_settings!(x.inner, verbose=Silent)
        end
    end

    if TimeLimit != nothing
        TimeLimit::Float
        x.settings[:time_limit] = !TimeLimit

        # Update time_limit setting if setup has already been performed
        if !x.perform_setup
            update_settings!(x.inner, time_limit=TimeLimit)
        end
    end



    x
end

# Do not update the problem instance settings using OSQP internal functions if the x is OSQPSolver and not OSQPMathProgModel
function MathProgBase.setparameters!(x::OSQPSolver; Silent = nothing, TimeLimit = nothing)
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent
    end
    if TimeLimit != nothing
        TimeLimit::Float
        x.settings[:time_limit] = TimeLimit
    end

    x
end

"""
    setwarmstart!(model::OSQPMathProgModel, x)

Warm start the solver solution with the primal variable `x`.

NB. The `warm_start!` function supports setting also the dual variables but in MathProgBase only
the primal variable warm start is supported. Note that OSQP performs warm starting automatically
when parts of the problem data change and a new optimize! is called.
"""
function MathProgBase.setwarmstart!(model::OSQPMathProgModel, x)
    copyto!(model.x_ws, x)
    model.run_warmstart = true
end


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
function MathProgBase.loadproblem!(model::OSQPMathProgModel, A, l, u, c, lb, ub, sense)

    (m, n) = size(A)  # Get problem dimensions
    initialize_dimensions!(model, n, m)  # Initialize problem dimensions

    if sense == :Min
        copyto!(model.q, c)
    elseif sense == :Max
        model.q .= .-c
    else
        error("Objective sense not recognized")
    end

    # Fix model sense
    model.sense = sense

    # Copy variables
    copyto!(model.lb, lb)
    copyto!(model.ub, ub)
    copyto!(model.l, l)
    copyto!(model.u, u)
    copyto!(model.A, sparse(A))  # Sparsify matrix A
    model
end

MathProgBase.getvarLB(model::OSQPMathProgModel) = model.l
function MathProgBase.setvarLB!(model::OSQPMathProgModel, l)
    # Copy new variable lower bounds
    copyto!(model.l, l)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, l=[model.lb; model.l])
    end

    # Reset problem solution
    resetproblem(model)

    return model
end

MathProgBase.getvarUB(model::OSQPMathProgModel) = model.u
function MathProgBase.setvarUB!(model::OSQPMathProgModel, u)
    # Copy new variable upper bounds
    copyto!(model.u, u)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, u=[model.ub; model.u])
    end

    # Reset problem solution
    resetproblem(model)

    return model
end


MathProgBase.getconstrLB(model::OSQPMathProgModel) = model.lb
function MathProgBase.setconstrLB!(model::OSQPMathProgModel, lb)
    # Copy new constraints lower bounds
    copyto!(model.lb, lb)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, l=[model.lb; model.l])
    end

    # Reset problem solution
    resetproblem(model)

    return model
end


MathProgBase.getconstrUB(model::OSQPMathProgModel) = model.ub
function MathProgBase.setconstrUB!(model::OSQPMathProgModel, ub)
    # Copy new constraints upper bounds
    copyto!(model.ub, ub)

    # Update upper bounds in the model
    if !model.perform_setup
        update!(model.inner, u=[model.ub; model.u])
    end

    # Reset problem solution
    resetproblem(model)

    return model
end


function MathProgBase.setobj!(model::OSQPMathProgModel, c)
    # Copy cost
    copyto!(model.q, c)

    # Negate cost if necessary
    if model.sense == :Max
        model.q .= .- model.q
    end

    if !model.perform_setup
        update!(model.inner, q=model.q)
    end

    # Reset problem solution
    resetproblem(model)

    return model

end

function MathProgBase.getobj(model::OSQPMathProgModel)
    if model.sense == :Min
        return model.q
    else
        return -model.q
    end
end


MathProgBase.getconstrmatrix(model::OSQPMathProgModel) = model.A


function MathProgBase.addvar!(model::OSQPMathProgModel, constridx, constrcoef, l, u, objcoef)

    # Get bounds if they are not set
    ((l == nothing) && (li = -Inf)) || (li = l)
    ((u == nothing) && (ui = Inf)) || (ui = u)

    # Change cost P, q
    model.q = [model.q; objcoef]
    if model.sense == :Max
        model.q[end] = -model.q[end]
    end
    model.P = blockdiag(model.P, spzeros(1,1))

    # Change constraints
    if !isempty(constrcoef) && !isempty(constridx)
        # Update A
        m = MathProgBase.numconstr(model)
        constr_xi = sparsevec(constridx, constrcoef, m)
        model.A = [model.A constr_xi]
    end

    # Update l and u
    model.l = [model.l; li]
    model.u = [model.u; ui]

    # Change ws
    model.x_ws = [model.x_ws; 0]

    # With new variable we need to perform setup
    model.perform_setup = true

    # Reset problem solution
    resetproblem(model)

end
MathProgBase.addvar!(model::OSQPMathProgModel, l, u, objcoef) = MathProgBase.addvar!(model, [], [], l, u, objcoef)

function MathProgBase.addconstr!(model::OSQPMathProgModel, varidx, coef, lb, ub)
    # Construct sparse vector with the new constraint
    n = MathProgBase.numvar(model)
    ai = sparsevec(varidx, coef, n)

    # Augment A, lb, ub
    if isempty(model.A)
        model.A = sparse(ai')
    else
        model.A = [model.A; ai']
    end
    model.lb = [model.lb; lb]
    model.ub = [model.ub; ub]

    # With the new constraint we must perform setup
    model.perform_setup = true

    # Reset problem solution
    resetproblem(model)

end

# TODO: delvars!
# TODO: delconstrs!
# TODO: changecoeffs!

MathProgBase.numlinconstr(model::OSQPMathProgModel) = MathProgBase.numconstr(model)
MathProgBase.getconstrsolution(model::OSQPMathProgModel) = model.A * MathProgBase.getsolution(model)
function MathProgBase.getreducedcosts(model::OSQPMathProgModel)
    checksolved(model)
    if model.sense == :Min
        return -model.results.y[end-(MathProgBase.numvar(model)-1):end]
    else
        return model.results.y[end-(MathProgBase.numvar(model)-1):end]
    end
end

function MathProgBase.getconstrduals(model::OSQPMathProgModel)
    checksolved(model)
    if model.sense == :Min
        return -model.results.y[1:MathProgBase.numconstr(model)]
    else
        return model.results.y[1:MathProgBase.numconstr(model)]
    end
end


function MathProgBase.getinfeasibilityray(model::OSQPMathProgModel)
    checksolved(model)

    if model.results.info.status in [:Primal_infeasible, :Primal_infeasible_inaccurate]
        m = MathProgBase.numconstr(model)
        # Returns infeasibility ray taking into account both bounds and constraints
        ray = -model.results.prim_inf_cert[1:MathProgBase.numconstr(model)]
    else
        error("Problem not infeasible")
    end
end


function MathProgBase.getunboundedray(model::OSQPMathProgModel)
    checksolved(model)

    if model.results.info.status in [:Dual_infeasible, :Dual_infeasible_inaccurate]
        return model.results.dual_inf_cert
    else
        error("Problem not unbounded")
    end
end

# TODO: getbasis

# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
MathProgBase.numquadconstr(model::OSQPMathProgModel) = 0

# setquadobj!(model::OSQPMathProgModel, Q::Matrix) = ((Qi, Qj, Qx) = findnz(Q); setquadobj!(model, Qi, Qj, Qx))
function MathProgBase.setquadobj!(model::OSQPMathProgModel, rowidx::Vector, colidx::Vector, quadval::Vector)
    nterms = length(quadval)
    @boundscheck length(rowidx) == nterms || error()
    @boundscheck length(colidx) == nterms || error()

    # Check if only the values have changed
    Pi, Pj, Px = findnz(model.P)
    if (rowidx == Pi) & (colidx == Pj) & !model.perform_setup
        if model.sense == :Max
            # Update matrix in MathProgBase model
            copyto!(model.P.nzval, -Px)
        else
            # Update matrix in MathProgBase model
            copyto!(model.P.nzval, Px)
        end
        # Update only nonzeros of P
        update!(model.inner, Px=model.P.nzval)

        # Reset solution status
        resetproblem(model)

        return model
    end

    # Create sparse matrix from indices
    # zero out coeffs
    for i = 1 : nterms
        @inbounds row = rowidx[i]
        @inbounds col = colidx[i]
        model.P[row, col] = 0
    end

    # add new coeffs and symmetrize
    for i = 1 : nterms
        @inbounds row = rowidx[i]
        @inbounds col = colidx[i]
        @inbounds val = quadval[i]
        model.P[row, col] += val
        model.P[col, row] = model.P[row, col]
    end

    # Change sign if maximizing
    if model.sense == :Max
        model.P .= .-model.P
    end

    # Need new setup when we set a new P
    model.perform_setup = true

    # Reset problem solution
    resetproblem(model)

    model
end

# Note: skipping quadconstr methods

end # module
