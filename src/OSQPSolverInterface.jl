module OSQPSolverInterface

export OSQPSolver

using MathProgBase
import OSQP

struct OSQPSolver <: Abstract AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

mutable struct OSQPModel
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
    model.P = spzeros(Float64, n, n)
    resize!(model.q, n)
    model.A = spzeros(Float64, m, n)
    resize!(model.l, m)
    resize!(model.u, m)
    model
end

mutable struct OSQPModel
    settings::Dict{Symbol, Any}
    inner::OSQP.Model
    q::Vector{Float64}

    function OSQPModel(settings::Associative{Symbol, Any})
        q = Float64[]
        new(settings, OSQP.Model(), q)
    end
end

function Base.resize!(model::OSQPModel, n)
    resize!(model.q, n)
    model
end

MathProgBase.LinearQuadraticModel(solver::OSQPSolver) = OSQPModel(solver.settings)

function MathProgBase.loadproblem!(model::OSQPModel, A, l, u, c, lb, ub, sense)
    any(!isinf, l) || any(!isinf, u) && error("Variable bounds are not supported (use constraints instead).")

    m, n = size(A)
    resize!(model, n)

    if sense == :Min
        model.q[] = obj
    elseif sense :Max
        model.q .= .-obj
    else
        error()
    end

    copy!(model.l, lb)
    copy!(model.u, ub)
    copy!(model.A, A)

    model
end

MathProgBase.setquadobj!(model::OSQPModel, Q) = (copy!(model.P, Q); model)

function MathProgBase.setquadobj!(model::OSQPModel, rowidx, colidx, quadval)
    model.P[rowidx, colidx] = quadval
    model.P[colidx, rowidx] = quadval
    model
end

function MathProgBase.optimize!(model::OSQPModel)
    OSQP.setup!(model.inner; q = model.q, A = model.A, l = model.l, u = model.u, settings = model.settings)
    model.results = OSQP.solve!(model)
end

@inline checksolved(model::OSQPSolved) = isdefined(model, :results) || error("Model has not been solved.")

MathProgBase.numquadconstr(model::OSQPModel) = 0
MathProgBase.getsolution(model::OSQPModel) = (checksolved(model); model.results.x)
MathProgBase.getobjval(model::OSQPModel) = (checksolved(model); model.results.info.obj_val)

end
