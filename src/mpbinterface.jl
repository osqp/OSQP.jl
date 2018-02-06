module OSQPMathProgBaseInterface

import MathProgBase
import OSQP

import MathProgBase: numvar, numconstr

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

function get_qp_dimensions(P::SparseMatrixCSC, q::Vector{Float64}, A::SparseMatrixCSC)
    # Check problem dimensions
    if P == nothing
        if q != nothing
            n = length(q)
        elseif A != nothing
            n = size(A, 2)
        end

    else
        n = size(P, 1)
    end

    if A == nothing
        m = 0
    else
        m = size(A, 1)
    end

    return n, m
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

# TODO: Remove results when the problem is updated 
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
# TODO: Fix sense!
function setsense!(model::OSQPMathProgModel, sense::Symbol) 
    if sense in [:Min, :Max]
        model.sense = sense
    else
        error("sense not recognized")
    end
end
getsense(model::OSQPMathProgModel) = model.sense

# TODO: Not necessarily right the following two lines
numvar(model::OSQPMathProgModel) = size(model.A, 2)
numconstr(model::OSQPMathProgModel) = size(model.A, 1)
freemodel!(model::OSQPMathProgModel) = OSQP.clean!(model.inner)

# NB COpy not implemented
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

# setparameters! unused at the moment. They are loaded using the initial loadproblem! call
# function setparameters!(x::Union{OSQPSolver, OSQPMathProgModel}; Silent = nothing)
#     if Silent != nothing
#         Silent::Bool
#         x.settings[:verbose] = !Silent
#     end
#     x
# end

# Check SCS.jl function to set also dual variables
setwarmstart!(model::OSQPMathProgModel, v) = OSQP.warm_start!(model.inner, x = v)


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
# TODO: loadproblem!(m::AbstractLinearQuadraticModel, filename::String)
function loadproblem!(model::OSQPMathProgModel, A, l, u, c, lb, ub, sense)
    (any(x -> x != -Inf, l) || any(x -> x != Inf, u)) && variablebounderror()
    if sense == :Min
        copy!(model.q, c)
    elseif sense == :Max
        model.q .= .-c
    else
        error("Objective sense not recognized")
    end

    m, n = size(A)
    resize!(model, n, m)
    copy!(model.l, lb)
    copy!(model.u, ub)
    copy!(model.A, A)
    model
end

# TODO: writeproblem
# TODO: Fix variable bounds
getvarLB(model::OSQPMathProgModel) = fill(-Inf, numvar(model))
setvarLB!(model::OSQPMathProgModel, l) = variablebounderror()
getvarUB(model::OSQPMathProgModel) = fill(Inf, numvar(model))
setvarUB!(model::OSQPMathProgModel, l) = variablebounderror()
getconstrLB(model::OSQPMathProgModel) = model.l
setconstrLB!(model::OSQPMathProgModel, lb) = (copy!(model.l, lb); model)
getconstrUB(model::OSQPMathProgModel) = model.u
setconstrUB!(model::OSQPMathProgModel, ub) = (copy!(model.u, ub); model)
setobj!(model::OSQPMathProgModel, c) = (copy!(model.q, c); model)
getconstrmatrix(model::OSQPMathProgModel) = model.A
# TODO: addvar!
# TODO: delvars!
# TODO: addconstr!
# TODO: delconstrs!
# TODO: changecoeffs!
numlinconstr(model::OSQPMathProgModel) = numconstr(model)
getconstrsolution(model::OSQPMathProgModel) = model.A * getsolution(model)
# TODO: getreducedcosts
getconstrduals(model::OSQPMathProgModel) = (checksolved(model); model.results.y)
# TODO: getinfeasibilityray
# TODO: getbasis
# TODO: getunboundedray


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
numquadconstr(model::OSQPMathProgModel) = 0

# TODO Change setquadobj! with proper sense
setquadobj!(model::OSQPMathProgModel, Q::Matrix) = (copy!(model.P, Q); model)
function setquadobj!(model::OSQPMathProgModel, rowidx::Vector, colidx::Vector, quadval::Vector)
    nterms = length(quadval)
    @boundscheck length(rowidx) == nterms || error()
    @boundscheck length(colidx) == nterms || error()

    # zero out coeffs
    for i = 1 : nterms
        @inbounds row = rowidx[i]
        @inbounds col = colidx[i]
        model.P[row, col] = 0
    end

    # add new coeffs
    for i = 1 : nterms
        @inbounds row = rowidx[i]
        @inbounds col = colidx[i]
        @inbounds val = quadval[i]
        model.P[row, col] += val
        model.P[col, row] = model.P[row, col]
    end

    model
end

# Note: skipping quadconstr methods

end # module
