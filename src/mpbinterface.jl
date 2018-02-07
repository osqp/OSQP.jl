module OSQPMathProgBaseInterface

import MathProgBase
import OSQP

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
        new(settings, Model(), perform_setup, :Min, P, q, A, lb, ub, l, u) # leave other fields undefined
    end
end

# TODO: Change here nothing to something more meaningful on the dimensions
function get_qp_variables(model::OSQPMathProgModel)
    # Check problem dimensions
    if model.P == nothing
        if model.q != nothing
            n = length(model.q)
        elseif model.A != nothing
            n = size(model.A, 2)
        end

    else
        n = size(model.P, 1)
    end

    return n
end

function get_qp_constraints(model::OSQPMathProgModel)
    if model.A == nothing
        m = 0
    else
        m = size(model.A, 1)
    end

    return m
end



# TODO: Do we need this?
# function Base.resize!(model::OSQPMathProgModel, n, m)
#     nchange = n != numvar(model)
#     mchange = m != numconstr(model)
#     if nchange
#         model.P = spzeros(Float64, n, n)
#         resize!(model.q, n)
#     end
#     if mchange
#         resize!(model.l, m)
#         resize!(model.u, m)
#     end
#     if nchange || mchange
#         model.A = spzeros(Float64, m, n)
#     end

#     model
# end

# Maybe add resetproblem function to remove results and reset perform_setup?
checksolved(model::OSQPMathProgModel) = isdefined(model, :results) || error("Model has not been solved.")

# Reset problem when setup has to be performed again
resetproblem(model::OSQPMathProgModel) = (model.results = nothing; model.perform_setup=true)

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
LinearQuadraticModel(solver::OSQPSolver) = OSQPMathProgModel(solver.settings)
getsolution(model::OSQPMathProgModel) = (checksolved(model); model.results.x)
getobjval(model::OSQPMathProgModel) = (checksolved(model); model.results.info.obj_val)

function optimize!(model::OSQPMathProgModel)
    # Perform setup only if necessary
    if model.perform_setup
        (n, _) = get_qp_dimensions(model.P, model.q, model.A)
        setup!(model.inner; P = model.P, q = model.q, 
               A = [model.A; sparse(I, n)], 
               l = [model.lb; model.l], 
               u = [model.ub; model.u], 
               settings = model.settings)
        model.perform_setup = false  # No longer needed to perform setup
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
    osqpstatus == :Dual_infeasible && (status = :DualityFailure)
    osqpstatus == :Dual_infeasible_inaccurate && (status = :DualityFailure)

    return status
end

# TODO: getobjbound
# TODO: getobjgap
getrawsolver(model::OSQPMathProgModel) = model.inner
getsolvetime(model::OSQPMathProgModel) = (checksolved(model); model.results.info.run_time)


function setsense!(model::OSQPMathProgModel, sense::Symbol) 
    if sense in [:Min, :Max]
        if sense == :Max & model.sense != :Max
            model.q *= -1
            model.P *= -1

            # Update problem data
            if !problem.perform_setup
                update!(model.inner, q=model.q)
                update!(model.inner, Px=triu(model.P.nzval))
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
freemodel!(model::OSQPMathProgModel) = OSQP.clean!(model.inner)

# NB Copy not implemented
# function Base.copy(model::OSQPMathProgModel)
#     ret = OSQPMathProgModel(model.settings)
#     resize!(ret, numvar(model), numconstr(model))
#     copy!(ret.P, model.P)
#     copy!(ret.q, model.q)
#     copy!(ret.A, model.A)
#     copy!(ret.l, model.l)
#     copy!(ret.u, model.u)
#     ret
# end

setvartype!(model::OSQPMathProgModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
getvartype(model::OSQPMathProgModel) = fill(:Cont, numvar(model))

#TODO: Fix setparameters!
# setparameters! unused at the moment. They are loaded using the initial loadproblem! call
# function setparameters!(x::Union{OSQPSolver, OSQPMathProgModel}; Silent = nothing)
#     if Silent != nothing
#         Silent::Bool
#         x.settings[:verbose] = !Silent
#     end
#     x
# end

# TODO: Check SCS.jl function to set also dual variables
setwarmstart!(model::OSQPMathProgModel, v) = OSQP.warm_start!(model.inner, x = v)


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
function loadproblem!(model::OSQPMathProgModel, A, l, u, c, lb, ub, sense)
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
    model.sense == :Max && model.q = -model.q   

    if !model.perform_setup
        update!(model.inner, q=model.q)
    end

    return model

end


getconstrmatrix(model::OSQPMathProgModel) = model.A
# TODO: addvar!
# TODO: delvars!
# TODO: addconstr!
# TODO: delconstrs!
# TODO: changecoeffs!
numlinconstr(model::OSQPMathProgModel) = numconstr(model)
getconstrsolution(model::OSQPMathProgModel) = model.A * getsolution(model)
getreducedcosts(model::OSQPMathProgModel) = (checksolved(model); model.results.y[end-get_qp_variables(model):end])
getconstrduals(model::OSQPMathProgModel) = (checksolved(model); model.results.y[1:get_qp_variables(model)])


function getinfeasibilityray(model::OSQPMathProgModel)
    checksolved(model)
    
    if model.info.status in [:Primal_infeasible, :Primal_infeasible_inaccurate]
        # NB This returns an infeasibility ray taking into account also the variable bounds
        return model.results.prim_inf_cert
    else
        error("Problem not infeasible")
    end
end


function getunboundedray(model::OSQPMathProgModel)
    checksolved(model)
    
    if model.info.status in [:Dual_infeasible, :Dual_infeasible_inaccurate]
        return model.results.dual_inf_cert
    else
        error("Problem not unbounded")
    end
end

# TODO: getbasis

# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
numquadconstr(model::OSQPMathProgModel) = 0

setquadobj!(model::OSQPMathProgModel, Q::Matrix) = ((Qi, Qj, Qx) = findnz(Q); setquadobj!(model, Qi, Qj, Qx))
function setquadobj!(model::OSQPMathProgModel, rowidx::Vector, colidx::Vector, quadval::Vector)
    nterms = length(quadval)
    @boundscheck length(rowidx) == nterms || error()
    @boundscheck length(colidx) == nterms || error()

    # Check if only the values have changed
    if isdefined(model, :P)
        Pi, Pj, Px = findnz(model.P)
        if rowidx == Pi & colidx == Pj & !model.perform_setup
            if model.sense == :Max
                # Update only nonzeros of P
                update!(model.inner, Px=-Px)
            else
                # Update only nonzeros of P
                update!(model.inner, Px=Px)         
            end
            return model       
        end
    end
    
    # Create sparse matrix from indices
    model.P = sparse(rowidx, colidx, quadval)

    if model.sense == :Max
        model.P = -model.P
    end
    
    # Need new setup when we set a new P
    model.perform_setup = true

    model
end

# Note: skipping quadconstr methods

end # module
