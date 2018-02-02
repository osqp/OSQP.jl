# Standard MathProg interrface
importall MathProgBase.SolverInterface

export OSQPSolver


struct OSQPSolver <: AbstractMathProgSolver
    settings::Dict{Symbol,Any}
end

OSQPSolver() = OSQPSolver(Dict{Symbol, Any}())

mutable struct OSQPMathProgModel <: AbstractLinearQuadraticModel
    # Problem settings
    settings::Dict{Symbol, Any}

    # Inner OSQP model structure
    inner::Model

    # Setup needed
    perform_setup::Bool

    # Problem data
    P::SparseMatrixCSC{Float64, Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int}
    l::Vector{Float64}
    u::Vector{Float64}

    # TODO: Write data in linprog format so that variable bounds are supported

    # Results
    results::Results

    function OSQPMathProgModel(settings::Associative{Symbol, Any})
        P = spzeros(Float64, 0, 0)
        q = Float64[]
        A = spzeros(Float64, 0, 0)
        l = Float64[]
        u = Float64[]
        perform_setup = true
        new(settings, Model(), perform_setup, P, q, A, l, u) # leave other fields undefined
    end
end

function Base.resize!(model::OSQPMathProgModel, n, m)
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

# TODO: Remove results when the problem is updated 
checksolved(model::OSQPMathProgModel) = isdefined(model, :results) || error("Model has not been solved.")
variablebounderror() = error("Variable bounds are not supported (use constraints instead).")
resetproblem(model::OSQPMathProgModel) = 

# http://mathprogbasejl.readthedocs.io/en/latest/solverinterface.html
LinearQuadraticModel(solver::OSQPSolver) = OSQPMathProgModel(solver.settings)
getsolution(model::OSQPMathProgModel) = (checksolved(model); model.results.x)
getobjval(model::OSQPMathProgModel) = (checksolved(model); model.results.info.obj_val)

function optimize!(model::OSQPMathProgModel)
    # Perform setup only if necessary
    if model.perform_setup
        setup!(model.inner; P = model.P, q = model.q, A = model.A, l = model.l, u = model.u, settings = model.settings)
        model.perform_setup = false  # No longer needed to perform setup
    end

    # Solve the problem and store the results
    model.results = solve!(model.inner)
end

function status(model::OSQPMathProgModel)::Symbol
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
    status == :Dual_infeasible_inaccurate && (ret = :DualityFailure)
    return ret
end

# TODO: getobjbound
# TODO: getobjgap
getrawsolver(model::OSQPMathProgModel) = model.inner
getsolvetime(model::OSQPMathProgModel) = (checksolved(model); model.results.info.run_time)
# TODO: setsense!
# TODO: getsense
numvar(model::OSQPMathProgModel) = size(model.A, 2)
numconstr(model::OSQPMathProgModel) = size(model.A, 1)
# TODO: freemodel!
# TODO: copy
setvartype!(model::OSQPMathProgModel, v::Vector{Symbol}) = any(x -> x != :Cont, v) && error("OSQP only supports continuous variables.")
getvartype(model::OSQPMathProgModel) = fill(:Cont, numvar(model))

# Removed setparameters! since it is unused at the moment
setwarmstart!(model::OSQPMathProgModel, v) = warm_start!(model.inner, x = v)


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#linearquadratic-models
# TODO: loadproblem!(m::AbstractLinearQuadraticModel, filename::String)
function loadproblem!(model::OSQPMathProgModel, A, l, u, c, lb, ub, sense)
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

# TODO: writeproblem(m::AbstractLinearQuadraticModel, filename::String)

# Edit problem data
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
# TODO: getsimplexiter
# TODO: getbarrieriter


# http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#quadratic-programming
numquadconstr(model::OSQPMathProgModel) = 0
setquadobj!(model::OSQPMathProgModel, Q) = (copy!(model.P, Q); model)

function setquadobj!(model::OSQPMathProgModel, rowidx, colidx, quadval)
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