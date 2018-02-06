module OSQPMathProgBaseInterface

import MathProgBase
import OSQP

import MathProgBase: numvar, numconstr

struct OSQPSolver <: MathProgBase.AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

OSQPSolver(; kwargs...) = OSQPSolver(Dict{Symbol, Any}(k => v for (k, v) in kwargs))

mutable struct OSQPModel <: MathProgBase.AbstractLinearQuadraticModel
    settings::Dict{Symbol, Any}
    inner::OSQP.Model

    P::SparseMatrixCSC{Float64, Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    l::Vector{Float64}
    u::Vector{Float64}

    xwarmstart::Vector{Float64}
    dowarmstart::Bool

    results::OSQP.Results

    function OSQPModel(settings::Associative{Symbol, Any})
        P = spzeros(Float64, 0, 0)
        q = Float64[]
        A = spzeros(Float64, 0, 0)
        l = Float64[]
        u = Float64[]
        xwarmstart = Float64[]
        dowarmstart = false
        new(copy(settings), OSQP.Model(), P, q, A, l, u, xwarmstart, dowarmstart) # leave results undefined
    end
end

function Base.resize!(model::OSQPModel, n, m)
    nchange = n != numvar(model)
    mchange = m != numconstr(model)
    if nchange
        model.P = spzeros(Float64, n, n)
        resize!(model.q, n)
        resize!(model.xwarmstart, n)
        model.dowarmstart = false
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
senseerror() = error("Only objective sense :Min is currently supported to avoid confusing behavior.")

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
MathProgBase.LinearQuadraticModel(solver::OSQPSolver) = OSQPModel(solver.settings)
MathProgBase.getsolution(model::OSQPModel) = (checksolved(model); model.results.x)
MathProgBase.getobjval(model::OSQPModel) = (checksolved(model); model.results.info.obj_val)

function MathProgBase.optimize!(model::OSQPModel)
    OSQP.setup!(model.inner; P = model.P, q = model.q, A = model.A, l = model.l, u = model.u, model.settings...)
    model.dowarmstart && OSQP.warm_start!(model.inner, x = model.xwarmstart)
    model.results = OSQP.solve!(model.inner)
end

function MathProgBase.status(model::OSQPModel)::Symbol
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

    return status
end

# TODO: getobjbound
# TODO: getobjgap
MathProgBase.getrawsolver(model::OSQPModel) = model.inner
MathProgBase.getsolvetime(model::OSQPModel) = (checksolved(model); model.results.info.run_time)
MathProgBase.setsense!(model::OSQPModel, sense::Symbol) = sense == :Min || senseerror()
MathProgBase.getsense(model::OSQPModel) = :Min
MathProgBase.numvar(model::OSQPModel) = size(model.A, 2)
MathProgBase.numconstr(model::OSQPModel) = size(model.A, 1)
MathProgBase.freemodel!(model::OSQPModel) = OSQP.clean!(model.inner)

function Base.copy(model::OSQPModel)
    ret = OSQPModel(model.settings)
    resize!(ret, numvar(model), numconstr(model))
    copy!(ret.P, model.P)
    copy!(ret.q, model.q)
    copy!(ret.A, model.A)
    copy!(ret.l, model.l)
    copy!(ret.u, model.u)
    ret
end

MathProgBase.setvartype!(model::OSQPModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
MathProgBase.getvartype(model::OSQPModel) = fill(:Cont, numvar(model))

function MathProgBase.setparameters!(x::Union{OSQPSolver, OSQPModel}; Silent = nothing)
    if Silent != nothing
        Silent::Bool
        x.settings[:verbose] = !Silent
    end
    x
end

MathProgBase.setwarmstart!(model::OSQPModel, v) = (copy!(model.xwarmstart, v); model.dowarmstart = true)


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
# TODO: loadproblem!(m::AbstractLinearQuadraticModel, filename::String)

function MathProgBase.loadproblem!(model::OSQPModel, A, l, u, c, lb, ub, sense)
    (any(x -> x != -Inf, l) || any(x -> x != Inf, u)) && variablebounderror()
    sense == :Min || senseerror()
    m, n = size(A)
    resize!(model, n, m)
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


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
MathProgBase.numquadconstr(model::OSQPModel) = 0
MathProgBase.setquadobj!(model::OSQPModel, Q::Matrix) = (copy!(model.P, Q); model)

function MathProgBase.setquadobj!(model::OSQPModel, rowidx::Vector, colidx::Vector, quadval::Vector)
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
