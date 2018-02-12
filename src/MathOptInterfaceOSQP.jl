module MathOptInterfaceOSQP

export OSQPOptimizer

using Compat
using MathOptInterface
using MathOptInterfaceUtilities

const MOI = MathOptInterface
const MOIU = MathOptInterfaceUtilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

import OSQP
import MathOptInterfaceUtilities: IndexMap

mutable struct OSQPOptimizer <: MOI.AbstractOptimizer # TODO: name?
    inner::OSQP.Model
    results::Union{OSQP.Results, Nothing}
    isempty::Bool
    # TODO: constant term?

    OSQPOptimizer() = new(OSQP.Model(), nothing, true)
end

hasresults(optimizer::OSQPOptimizer) = optimizer.results != nothing

function MOI.empty!(optimizer::OSQPOptimizer)
    optimizer.inner = OSQP.Model()
    optimizer.results = nothing
    optimizer.isempty = true
end

MOI.isempty(optimizer::OSQPOptimizer) = optimizer.isempty

function MOI.copy!(dest::OSQPOptimizer, src::MOI.ModelLike)
    # TODO: cangets (awaiting https://github.com/JuliaOpt/MathOptInterface.jl/pull/210)

    # Type aliases
    Affine = MOI.ScalarAffineFunction{Float64}
    SingleVariable = MOI.SingleVariable
    Quadratic = MOI.ScalarQuadraticFunction{Float64}
    Interval = MOI.Interval{Float64}

    # Empty
    MOI.empty!(dest)

    # Check constraint types
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        MOI.supportsconstraint(dest, F, S) || return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Unsupported $F-in-$S constraint", IndexMap())
    end

    # Set up index map
    idxmap = IndexMap() # TODO: could even use a typed version of this, as there's only one type of constraint right now.
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = VI(i)
    end
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{Affine, Interval}()) # TODO: canget
    for i in eachindex(cis_src)
        idxmap[cis_src[i]] = CI{Affine, Interval}(i)
    end

    # Get sizes
    n = MOI.get(src, MOI.NumberOfVariables())
    m = MOI.get(src, MOI.NumberOfConstraints{Affine, Interval}())

    # Allocate storage for problem data
    P = spzeros(n, n)
    q = zeros(n)
    A = spzeros(m, n)
    l = zeros(m)
    u = zeros(m)

    # Process objective function
    # prefer the simplest form
    if MOI.canget(src, MOI.ObjectiveFunction{SingleVariable}())
        processobjective!(P, q, MOI.get(src, MOI.ObjectiveFunction{SingleVariable}()), idxmap)
    elseif MOI.canget(src, MOI.ObjectiveFunction{Affine}())
        processobjective!(P, q, MOI.get(src, MOI.ObjectiveFunction{Affine}()), idxmap)
    elseif MOI.canget(src, MOI.ObjectiveFunction{Quadratic}())
        processobjective!(P, q, MOI.get(src, MOI.ObjectiveFunction{Quadratic}()), idxmap)
    else
        return MOI.CopyResult(MOI.CopyOtherError, "No suitable objective function found", idxmap)
    end
    sense = MOI.get(src, MOI.ObjectiveSense())
    @assert sense == MOI.MinSense || sense == MOI.FeasibilitySense
    # TODO: constant term

    # Process constraints
    processconstraints!(A, l, u, src, idxmap, Affine, Interval)

    # Load data into OSQP Model
    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u) # TODO: settings

    # Process optimizer attributes
    # TODO

    # Process variable attributes
    if MOI.canget(src, MOI.VariablePrimalStart(), VI)
        x = zeros(n)
        i = 1
        for vi in vis_src
            x[i] = get(src, MOI.VariablePrimalStart(), vi)
            i += 1
        end
        OSQP.warm_start!(src, x = x)
    end

    if MOI.canget(src, MOI.ConstraintDualStart(), CI{Affine, Interval})
        y = zeros(m)
        i = 1
        for ci in cis_src
            y[i] = get(src, MOI.ConstraintDualStart(), ci)
            i += 1
        end
    end

    # Finish up
    dest.isempty = false
    return MOI.CopyResult(MOI.CopySuccess, "", idxmap)
end

# These could probably move to MOIU:
function process_affine_objective!(P::SparseMatrixCSC, q::Vector, variables::Vector{VI}, coefficients::Vector, idxmap)
    q[:] = 0
    n = length(coefficients)
    @assert length(variables) == n
    for i = 1 : n
        q[idxmap[var[i]].value] += coefficients[i]
    end
end

function processobjective!(P::SparseMatrixCSC, q::Vector, objfun::MOI.ScalarAffineFunction, idxmap)
    process_affine_objective!(P, q, objfun.variables, objfun.coefficients, idxmap)
end

function processobjective!(P::SparseMatrixCSC, q::Vector, objfun::MOI.ScalarQuadraticFunction, idxmap)
    process_affine_objective!(P, q, objfun.affine_variables, objfun.affine_coefficients, idxmap)
    nquadratic = length(objfun.quadratic_coefficients)
    @assert length(objfun.quadratic_rowvariables) == length(objfun.quadratic_colvariables) == nquadratic
    P[:] = 0
    for i = 1 : nquadratic
        row = idxmap[objfun.quadratic_rowvariables[i]].value
        col = idxmap[objfun.quadratic_colvariables[i]].value
        P[row, col] += objfun.quadratic_coefficients[i]
        P[col, row] = P[row, col]
    end
end

function processobjective!(P::SparseMatrixCSC, q::Vector, objfun::MOI.SingleVariable, idxmap)
    q[:] = 0
    P[:] = 0
    q[idxmap[objfun.variable.value]] = 1
    nothing
end

function processconstraints!(A::SparseMatrixCSC, l::Vector, u::Vector, src::MOI.ModelLike, idxmap, F::Type{<:MOI.ScalarAffineFunction}, S::Type{<:MOI.Interval})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    row = 1
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        l[row] = s.lower - f.constant
        u[row] = s.upper - f.constant
        n = length(f.coefficients)
        @assert length(f.variables) == n
        for i = 1 : n
            var = f.variables[i]
            coeff = f.coefficients[i]
            col = idxmap[var].value
            A[row, col] = coeff
        end
        row += 1
    end
end


## Optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = MOI.MinSense

MOI.canget(optimizer::OSQPOptimizer, ::MOI.NumberOfVariables) = !optimizer.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
MOI.get(optimizer::OSQPOptimizer, ::MOI.NumberOfVariables) = OSQP.dimensions(optimizer.model)[1]

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = MOI.get(optimizer, MOI.NumberOfVariables())
MOI.get(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = [VI(i) for i = 1 : get(optimizer, MOI.NumberOfVariables())] # TODO: support for UnitRange would be nice

MOI.canget(optimizer::OSQPOptimizer, ::MOI.NumberOfConstraints) = !optimizer.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
MOI.get(optimizer::OSQPOptimizer, ::MOI.NumberOfConstraints) = OSQP.dimensions(optimizer.model)[2]

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfConstraints) = false # TODO
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfConstraintIndices) = false # TODO
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfModelAttributesSet) = false # currently not exposed
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfVariableAttributesSet) = false # currently not exposed (for warmstart, rest is N/A)
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfConstraintAttributesSet) = false # currently not exposed

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintPrimal, ::CI) = false # constraint primal (or A matrix) not currently exposed by OSQP interface


## Solver optimizer:
MOI.optimize!(optimizer::OSQPOptimizer) = (optimizer.results = OSQP.solve!(optimizer.inner))
MOI.free!(optimizer::OSQPOptimizer) = OSQP.clean!(optimizer.inner)


## Solver optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.RawSolver) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.RawSolver) = optimizer.inner

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ResultCount) = hasresults(optimizer) # TODO: or true?
MOI.get(optimizer::OSQPOptimizer, ::MOI.ResultCount) = 1

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveFunction) = false # currently not exposed

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveValue) = hasresults(optimizer)
MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveValue) = optimizer.results.info.obj_val

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveBound) = false # currently not exposed
MOI.canget(optimizer::OSQPOptimizer, ::MOI.RelativeGap) = false # currently not exposed

MOI.canget(optimizer::OSQPOptimizer, ::MOI.SolveTime) = hasresults(optimizer)
MOI.get(optimizer::OSQPOptimizer, ::MOI.SolveTime) = optimizer.results.info.run_time

MOI.canget(optimizer::OSQPOptimizer, ::MOI.TerminationStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, ::MOI.TerminationStatus)
    # Note that the :Dual_infeasible and :Primal_infeasible are mapped to MOI.Success
    # because OSQP can return a proof of infeasibility. For the same reason,
    # :Primal_infeasible_inaccurate is mapped to MOI.AlmostSuccess

    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif osqpstatus == :Interrupted
        MOI.Interrupted
    elseif osqpstatus == :Dual_infeasible
        MOI.Success
    elseif osqpstatus == :Primal_infeasible
        MOI.Success
    elseif osqpstatus == :Max_iter_reached
        MOI.IterationLimit
    elseif osqpstatus == :Solved
        MOI.Success
    elseif osqpstatus == :Solved_inaccurate
        MOI.AlmostSuccess
    elseif osqpstatus == :Primal_infeasible_inaccurate
        MOI.AlmostSuccess
    end
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.PrimalStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, ::MOI.PrimalStatus)
    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif osqpstatus == :Primal_infeasible
        MOI.InfeasibilityCertificate
    elseif osqpstatus == :Solved
        MOI.FeasiblePoint
    elseif osqpstatus == :Primal_infeasible_inaccurate
        MOI.NearlyInfeasibilityCertificate
    else # :Interrupted, :Dual_infeasible, :Max_iter_reached, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.DualStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, ::MOI.DualStatus)
    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif osqpstatus == :Dual_infeasible
        MOI.InfeasibilityCertificate
    elseif osqpstatus == :Solved
        MOI.FeasiblePoint
    else # :Interrupted, :Primal_infeasible, :Max_iter_reached, :Primal_infeasible_inaccurate, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end

# TODO: solver-specific attributes


## Variables and constraints:
MOI.candelete(optimizer::OSQPOptimizer, index::MOI.Index) = false
MOI.isvalid(optimizer::OSQPOptimizer, vi::VI) = vi.value ∈ 1 : get(optimizer, MOI.NumberOfVariables())
MOI.canaddvariable(optimizer::OSQPOptimizer) = false


## Variable attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.VariablePrimalStart, ::Type{VI}) = false # currently not exposed, but could be
MOI.canset(optimizer::OSQPOptimizer, ::MOI.VariablePrimalStart, ::Type{VI}) = false # TODO: need selective way of updating primal start

MOI.canget(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, ::Type{VI}) = hasresults(optimizer) && optimizer.results.info.status ∈ OSQP.SOLUTION_PRESENT
MOI.get(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, vi::VI) = optimizer.results.x[vi.value]
MOI.get(optimizer::OSQPOptimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(optimizer, a, vi) # TODO: copied from SCS. Necessary?


## Constraints:
MOI.isvalid(optimizer::OSQPOptimizer, ci::CI) = ci.value ∈ 1 : get(optimizer, MOI.NumberOfConstraints())
MOI.canaddconstraint(optimizer::OSQPOptimizer, ::Type{F}, ::Type{S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet} = false

# TODO: can't modifyconstraint! with AbstractFunction because selective way to update A not exposed
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, ::Type{MOI.ScalarAffineFunction}) = false
# MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, func::MOI.ScalarAffineFunction)

# TODO: can't modifyconstraint! with AbstractSet because selective way to update lb and ub not exposed
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, ::Type{MOI.Interval}) = false
# MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, set::MOI.Interval)

# TODO: partial change with MultirowChange

MOI.supportsconstraint(optimizer::OSQPOptimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true


## Constraint attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintPrimalStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canset(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{VI}) = false # TODO: need selective way of updating primal start
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintPrimal, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDual, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintFunction, ::Type{<:CI}) = false # TODO
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintSet, ::Type{<:CI}) = false # TODO


# Objective modification
MOI.canmodifyobjective(optimizer::OSQPOptimizer, ::Type{MOI.ScalarCoefficientChange}) = false # TODO: selective way of updating objective coefficients not exposed

end # module
