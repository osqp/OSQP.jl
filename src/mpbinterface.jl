module OSQPMathProgBaseInterface

using OSQP: Model, Results, setup!, solve!, update!, clean!, update_settings!, warm_start!
importall MathProgBase.SolverInterface

struct OSQPSolver <: AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

OSQPSolver(; kwargs...) = OSQPSolver(Dict{Symbol, Any}(k => v for (k, v) in kwargs))

mutable struct OSQPMathProgModel <: AbstractLinearQuadraticModel
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

    # Results
    results::Results

    function OSQPMathProgModel(settings::Associative{Symbol, Any})
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
        new(settings, Model(), perform_setup, 
            :Min, P, q, A, lb, ub, l, u, 
            x_ws, run_warmstart) # leave other fields undefined
    end
end

function get_qp_variables(model::OSQPMathProgModel)
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

get_qp_constraints(model::OSQPMathProgModel) = size(model.A, 1)

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
checksolved(model::OSQPMathProgModel) = isdefined(model, :results) || error("Model has not been solved.")

# Reset problem when setup has to be performed again
resetproblem(model::OSQPMathProgModel) = (model.results = nothing; model.perform_setup=true)

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
LinearQuadraticModel(solver::OSQPSolver) = OSQPMathProgModel(solver.settings)
getsolution(model::OSQPMathProgModel) = (checksolved(model); model.results.x)
function getobjval(model::OSQPMathProgModel) 
    checksolved(model)
    
    if model.sense == :Min
        return model.results.info.obj_val
    else
        return -model.results.info.obj_val
    end
end

function optimize!(model::OSQPMathProgModel)
    # Perform setup only if necessary
    if model.perform_setup
        n = get_qp_variables(model)
        setup!(model.inner; P = model.P, q = model.q, 
               A = [model.A; sparse(I, n)], 
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
end

function status(model::OSQPMathProgModel)::Symbol
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

    return status
end

# TODO: getobjbound
# TODO: getobjgap
getrawsolver(model::OSQPMathProgModel) = model.inner
getsolvetime(model::OSQPMathProgModel) = (checksolved(model); model.results.info.run_time)


function setsense!(model::OSQPMathProgModel, sense::Symbol) 
    if sense in [:Min, :Max]
        if (sense == :Max) & (model.sense != :Max)
            model.q *= -1
            model.P *= -1

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
end

getsense(model::OSQPMathProgModel) = model.sense
numvar(model::OSQPMathProgModel) = get_qp_variables(model)
numconstr(model::OSQPMathProgModel) = get_qp_constraints(model)
freemodel!(model::OSQPMathProgModel) = clean!(model.inner)

function Base.copy(model::OSQPMathProgModel)
    ret = OSQPMathProgModel(model.settings)
    initialize_dimensions!(ret, numvar(model), numconstr(model))
    ret.sense = model.sense
    copy!(ret.P, model.P)
    copy!(ret.q, model.q)
    copy!(ret.A, model.A)
    copy!(ret.lb, model.lb)
    copy!(ret.ub, model.ub)
    copy!(ret.l, model.l)
    copy!(ret.u, model.u)
    ret
end

setvartype!(model::OSQPMathProgModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
getvartype(model::OSQPMathProgModel) = fill(:Cont, numvar(model))

# Set solver independent parameters: Silent is the only supported one for now
function setparameters!(x::OSQPMathProgModel; Silent = nothing)
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent

        # Update silent setting if setup has already been performed
        if !x.perform_setup
            update_settings!(x.inner, verbose=Silent)
        end
    end

    x
end

# Do not update the problem instance settings using OSQP internal functions if the x is OSQPSolver and not OSQPMathProgModel
function setparameters!(x::OSQPSolver; Silent = nothing)  
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent
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
function setwarmstart!(model::OSQPMathProgModel, x)
    copy!(model.x_ws, x)
    model.run_warmstart = true
end


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
function loadproblem!(model::OSQPMathProgModel, A, l, u, c, lb, ub, sense)

    (m, n) = size(A)  # Get problem dimensions
    initialize_dimensions!(model, n, m)  # Initialize problem dimensions

    if sense == :Min
        copy!(model.q, c)
    elseif sense == :Max
        model.q .= .-c
    else
        error("Objective sense not recognized")
    end

    # Copy variables
    copy!(model.q, c)
    copy!(model.lb, lb)
    copy!(model.ub, ub)
    copy!(model.l, l)
    copy!(model.u, u)
    copy!(model.A, sparse(A))  # Sparsify matrix A
    model
end

getvarLB(model::OSQPMathProgModel) = model.l
function setvarLB!(model::OSQPMathProgModel, l)
    # Copy new variable lower bounds
    copy!(model.l, l)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, l=[model.lb; model.l])
    end
    
    return model
end

getvarUB(model::OSQPMathProgModel) = model.u
function setvarUB!(model::OSQPMathProgModel, u)
    # Copy new variable upper bounds
    copy!(model.u, u)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, u=[model.ub; model.u])
    end
    
    return model
end


getconstrLB(model::OSQPMathProgModel) = model.lb
function setconstrLB!(model::OSQPMathProgModel, lb)
    # Copy new constraints lower bounds
    copy!(model.lb, lb)

    # Update lower bounds in the model
    if !model.perform_setup
        update!(model.inner, l=[model.lb; model.l])
    end
    
    return model
end


getconstrUB(model::OSQPMathProgModel) = model.ub
function setconstrUB!(model::OSQPMathProgModel, ub)
    # Copy new constraints upper bounds
    copy!(model.ub, ub)

    # Update upper bounds in the model
    if !model.perform_setup
        update!(model.inner, u=[model.ub; model.u])
    end
    
    return model
end


function setobj!(model::OSQPMathProgModel, c)
    # Copy cost
    copy!(model.q, c)
    
    # Negate cost if necessary
    if model.sense == :Max
        model.q *= -1
    end

    if !model.perform_setup
        update!(model.inner, q=model.q)
    end

    return model

end

function getobj(model::OSQPMathProgModel) 
    if model.sense == :Min
        return model.q
    else
        return -model.q
    end
end


getconstrmatrix(model::OSQPMathProgModel) = model.A


function addvar!(model::OSQPMathProgModel, constridx, constrcoef, l, u, objcoef)

    # Get bounds if they are not set
    ((l == nothing) && (li = -Inf)) || (li = l)
    ((u == nothing) && (ui = Inf)) || (ui = u)

    # Change cost P, q
    if model.sense == :Min
        qi = objcoef
    else
        qi = -objcoef
    end
    model.q = [model.q; qi]
    model.P = blkdiag(model.P, spzeros(1,1))

    # Change constraints
    if !isempty(constrcoef) && !isempty(constridx)
        # Update A
        m = numconstr(model)
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

end
addvar!(model::OSQPMathProgModel, l, u, objcoef) = addvar!(model, [], [], l, u, objcoef)

function addconstr!(model::OSQPMathProgModel, varidx, coef, lb, ub)
    # Construct sparse vector with the new constraint
    n = numvar(model)
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

end

# TODO: delvars!
# TODO: delconstrs!
# TODO: changecoeffs!

numlinconstr(model::OSQPMathProgModel) = numconstr(model)
getconstrsolution(model::OSQPMathProgModel) = model.A * getsolution(model)
function getreducedcosts(model::OSQPMathProgModel) 
    checksolved(model)
    if model.sense == :Min
        return -model.results.y[end-(get_qp_variables(model)-1):end]
    else
        return model.results.y[end-(get_qp_variables(model)-1):end]
    end
end

function getconstrduals(model::OSQPMathProgModel) 
    checksolved(model)
    if model.sense == :Min
        return -model.results.y[1:get_qp_constraints(model)]
    else
        return model.results.y[1:get_qp_constraints(model)]
    end
end


function getinfeasibilityray(model::OSQPMathProgModel)
    checksolved(model)
    
    if model.results.info.status in [:Primal_infeasible, :Primal_infeasible_inaccurate]
        m = numconstr(model)
        # Returns infeasibility ray taking into account both bounds and constraints
        ray = -model.results.prim_inf_cert
    else
        error("Problem not infeasible")
    end
end


function getunboundedray(model::OSQPMathProgModel)
    checksolved(model)
    
    if model.results.info.status in [:Dual_infeasible, :Dual_infeasible_inaccurate]
        return model.results.dual_inf_cert
    else
        error("Problem not unbounded")
    end
end

# TODO: getbasis

# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
numquadconstr(model::OSQPMathProgModel) = 0

# setquadobj!(model::OSQPMathProgModel, Q::Matrix) = ((Qi, Qj, Qx) = findnz(Q); setquadobj!(model, Qi, Qj, Qx))
function setquadobj!(model::OSQPMathProgModel, rowidx::Vector, colidx::Vector, quadval::Vector)
    nterms = length(quadval)
    @boundscheck length(rowidx) == nterms || error()
    @boundscheck length(colidx) == nterms || error()

    # Check if only the values have changed
    if isdefined(model, :P)
        Pi, Pj, Px = findnz(model.P)
        if (rowidx == Pi) & (colidx == Pj) & !model.perform_setup
            if model.sense == :Max
                # Update matrix in MathProgBase model
                copy!(model.P.nzval, -Px)
            else
                # Update matrix in MathProgBase model
                copy!(model.P.nzval, Px)      
            end
            # Update only nonzeros of P
            update!(model.inner, Px=model.P.nzval)
            return model       
        end
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
        model.P *= -1
    end
    
    # Need new setup when we set a new P
    model.perform_setup = true

    model
end

# Note: skipping quadconstr methods

end # module
