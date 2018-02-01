module OSQPSolverInterface

using MathProgBase
import OSQP

import MathProgBase: numvar, numconstr

struct OSQPSolver <: MathProgBase.AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

OSQPSolver() = OSQPSolver(Dict{Symbol, Any}())

mutable struct OSQPModel <: MathProgBase.AbstractLinearQuadraticModel
    settings::Dict{Symbol, Any}
    inner::OSQP.Model

    P::SparseMatrixCSC{Float64, Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    l::Vector{Float64}
    u::Vector{Float64}

    results::OSQP.Results

    function OSQPModel(settings::Associative{Symbol, Any})
        P = spzeros(Float64, 0, 0)
        q = Float64[]
        A = spzeros(Float64, 0, 0)
        l = Float64[]
        u = Float64[]
        new(settings, OSQP.Model(), P, q, A, l, u) # leave results undefined
    end
end

function Base.resize!(model::OSQPModel, n, m)
    nchange = n != numvar(model)
    mchange = m != numconstr(model)
    if nchange
        model.P = spzeros(Float64, n, n)
        resize!(model.q, n)
    end
    if mchange
        resize!(model.l, m)
        resize!(model.u, m)
    end
    if nchange || mchange
        model.A = spzeros(Float64, m, n)
    end

    model
end

checksolved(model::OSQPModel) = isdefined(model, :results) || error("Model has not been solved.")
variablebounderror() = error("Variable bounds are not supported (use constraints instead).")

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
MathProgBase.LinearQuadraticModel(solver::OSQPSolver) = OSQPModel(solver.settings)
MathProgBase.getsolution(model::OSQPModel) = (checksolved(model); model.results.x)
MathProgBase.getobjval(model::OSQPModel) = (checksolved(model); model.results.info.obj_val)

function MathProgBase.optimize!(model::OSQPModel)
    OSQP.setup!(model.inner; P = model.P, q = model.q, A = model.A, l = model.l, u = model.u, settings = model.settings)
    model.results = OSQP.solve!(model.inner)
end

function MathProgBase.status(model::OSQPModel)::Symbol
    checksolved(model)
    status = model.results.info.status
    ret = status # if OSQP status can't be mapped to a standard return value, just return as is

    # map to standard return status values:
    status == :Solved && (ret = :Optimal)
    status == :Max_iter_reached && (ret = :UserLimit)
    status == :Interrupted && (ret = :UserLimit) # following Gurobi.jl
    status == :Primal_infeasible && (ret = :Infeasible)
    status == :Primal_infeasible_inaccurate && (ret = :Infeasible)
    status == :Dual_infeasible && (ret = :DualityFailure)

    return ret
end

# TODO: getobjbound
# TODO: getobjgap
MathProgBase.getrawsolver(model::OSQPModel) = model.inner
MathProgBase.getsolvetime(model::OSQPModel) = (checksolved(model); model.results.info.run_time)
# TODO: setsense!
# TODO: getsense
MathProgBase.numvar(model::OSQPModel) = size(model.A, 2)
MathProgBase.numconstr(model::OSQPModel) = size(model.A, 1)
# TODO: freemodel!
# TODO: copy
MathProgBase.setvartype!(model::OSQPModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
MathProgBase.getvartype(model::OSQPModel) = fill(:Cont, numvar(model))

function MathProgBase.setparameters!(x::Union{OSQPSolver, OSQPModel}; Silent = nothing)
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent
    end
end

MathProgBase.setwarmstart!(model::OSQPModel, v) = OSQP.warm_start!(model.inner, x = v)


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
# TODO: loadproblem!(m::AbstractLinearQuadraticModel, filename::String)

function MathProgBase.loadproblem!(model::OSQPModel, A, l, u, c, lb, ub, sense)
    (any(x -> x != -Inf, l) || any(x -> x != Inf, u)) && variablebounderror()
    m, n = size(A)
    resize!(model, n, m)

    if sense == :Min
        copy!(model.q, c)
    elseif sense == :Max
        model.q .= .-c
    else
        error("Objective sense not recognized")
    end

    copy!(model.l, lb)
    copy!(model.u, ub)
    copy!(model.A, A)

    model
end

# TODO: writeproblem
MathProgBase.getvarLB(model::OSQPModel) = fill(-Inf, numvar(model))
MathProgBase.setvarLB!(model::OSQPModel, l) = variablebounderror()
MathProgBase.getvarUB(model::OSQPModel) = fill(Inf, numvar(model))
MathProgBase.setvarUB!(model::OSQPModel, l) = variablebounderror()
MathProgBase.getconstrLB(model::OSQPModel) = model.l
MathProgBase.setconstrLB!(model::OSQPModel, lb) = (copy!(model.l, lb); model)
MathProgBase.getconstrUB(model::OSQPModel) = model.u
MathProgBase.setconstrUB!(model::OSQPModel, ub) = (copy!(model.u, ub); model)
MathProgBase.setobj!(model::OSQPModel, c) = (copy!(model.q, c); model)
MathProgBase.getconstrmatrix(model::OSQPModel) = model.A
# TODO: addvar!
# TODO: delvars!
# TODO: addconstr!
# TODO: delconstrs!
# TODO: changecoeffs!
MathProgBase.numlinconstr(model::OSQPModel) = numconstr(model)
MathProgBase.getconstrsolution(model::OSQPModel) = model.A * getsolution(model)
# TODO: getreducedcosts
MathProgBase.getconstrduals(model::OSQPModel) = (checksolved(model); model.results.y)
# TODO: getinfeasibilityray
# TODO: getbasis
# TODO: getunboundedray
# TODO: getsimplexiter
# TODO: getbarrieriter


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
MathProgBase.numquadconstr(model::OSQPModel) = 0
MathProgBase.setquadobj!(model::OSQPModel, Q) = (copy!(model.P, Q); model)

function MathProgBase.setquadobj!(model::OSQPModel, rowidx, colidx, quadval)
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
