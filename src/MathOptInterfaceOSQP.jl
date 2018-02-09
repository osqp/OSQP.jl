module MathOptInterfaceOSQP

export OSQPInstance

using Compat
using MathOptInterface
using MathOptInterfaceUtilities

const MOI = MathOptInterface
const MOIU = MathOptInterfaceUtilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

import OSQP
import MathOptInterfaceUtilities: IndexMap

mutable struct OSQPInstance <: MOI.AbstractSolverInstance
    inner::OSQP.Model
    results::Union{OSQP.Results, Nothing}
    isempty::Bool
    # TODO: constant term?

    OSQPInstance() = new(OSQP.Model(), nothing, true)
end

hasresults(instance::OSQPInstance) = instance.results != nothing

function MOI.empty!(instance::OSQPInstance)
    instance.inner = OSQP.Model()
    instance.results = nothing
    instance.isempty = true
end

MOI.isempty(instance::OSQPInstance) = instance.isempty

function MOI.copy!(dest::OSQPInstance, src::MOI.AbstractInstance)
    MOI.empty!(dest)

    # Set up index map
    idxmap = MOI.IndexMap() # TODO: could even get a typed one
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = VI(i)
    end
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices())
    for i in eachindex(cis_src)
        ci_src = cis_src[i]
        if ci_src isa CI{MOI.ScalarAffineFunction, MOI.Interval}
            idxmap[cis_src[i]] = CI{MOI.ScalarAffineFunction, MOI.Interval}(i)
        else
            return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Unsupported constraint of type $(typeof(ci_src))", idxmap)
        end
    end

    # Get sizes
    @assert MOI.canget(src, MOI.NumberOfVariables())
    n = MOI.get(src, MOI.NumberOfVariables())
    @assert MOI.canget(src, MOI.NumberOfConstraints())
    m = MOI.get(src, MOI.NumberOfConstraints())

    # Allocate storage for problem data
    P = spzeros(n, n)
    q = zeros(n)
    A = spzeros(m, n)
    l = zeros(m)
    u = zeros(m)

    # Process objective function
    @assert MOI.canget(src, MOI.ObjectiveFunction())
    obj = MOI.get(src, MOI.ObjectiveFunction())
    processobjective!(P, q, obj, idxmap)
    # TODO: constant term

    # Process instance attributes
    @assert MOI.get(src, MOI.ObjectiveSense()) == MOI.get(dest, MOI.ObjectiveSense())

    # Process variable attributes
    # TODO: VariablePrimalStart

    # Process constraint attributes
    # TODO: ConstraintPrimalStart
    # TODO: ConstraintDualStart


    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u) # TODO: settings

    dest.isempty = false

    return MOI.CopyResult(MOI.CopySuccess, "", idxmap)
end

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


## Instance attributes:
MOI.canget(instance::OSQPInstance, ::MOI.ObjectiveSense) = true
MOI.get(instance::OSQPInstance, ::MOI.ObjectiveSense) = MOI.MinSense

MOI.canget(instance::OSQPInstance, ::MOI.NumberOfVariables) = !instance.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
MOI.get(instance::OSQPInstance, ::MOI.NumberOfVariables) = OSQP.dimensions(instance.model)[1]

MOI.canget(instance::OSQPInstance, ::MOI.ListOfVariableIndices) = MOI.get(instance, MOI.NumberOfVariables())
MOI.get(instance::OSQPInstance, ::MOI.ListOfVariableIndices) = [VI(i) for i = 1 : get(instance, MOI.NumberOfVariables())] # TODO: support for UnitRange would be nice

MOI.canget(instance::OSQPInstance, ::MOI.NumberOfConstraints) = !instance.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
MOI.get(instance::OSQPInstance, ::MOI.NumberOfConstraints) = OSQP.dimensions(instance.model)[2]

MOI.canget(instance::OSQPInstance, ::MOI.ListOfConstraints) = false # TODO
MOI.canget(instance::OSQPInstance, ::MOI.ListOfConstraintIndices) = false # TODO
MOI.canget(instance::OSQPInstance, ::MOI.ListOfInstanceAttributesSet) = false # currently not exposed
MOI.canget(instance::OSQPInstance, ::MOI.ListOfVariableAttributesSet) = false # currently not exposed (for warmstart, rest is N/A)
MOI.canget(instance::OSQPInstance, ::MOI.ListOfConstraintAttributesSet) = false # currently not exposed

MOI.canget(instance::OSQPInstance, ::MOI.ConstraintPrimal, ::CI) = false # constraint primal (or A matrix) not currently exposed by OSQP interface


## Solver instance:
MOI.optimize!(instance::OSQPInstance) = (instance.results = OSQP.solve!(instance.inner))
MOI.free!(instance::OSQPInstance) = OSQP.clean!(instance.inner)


## Solver instance attributes:
MOI.canget(instance::OSQPInstance, ::MOI.RawSolver) = true
MOI.get(instance::OSQPInstance, ::MOI.RawSolver) = instance.inner

MOI.canget(instance::OSQPInstance, ::MOI.ResultCount) = hasresults(instance) # TODO: or true?
MOI.get(instance::OSQPInstance, ::MOI.ResultCount) = 1

MOI.canget(instance::OSQPInstance, ::MOI.ObjectiveFunction) = false # currently not exposed

MOI.canget(instance::OSQPInstance, ::MOI.ObjectiveValue) = hasresults(instance)
MOI.get(instance::OSQPInstance, ::MOI.ObjectiveValue) = instance.results.info.obj_val

MOI.canget(instance::OSQPInstance, ::MOI.ObjectiveBound) = false # currently not exposed
MOI.canget(instance::OSQPInstance, ::MOI.RelativeGap) = false # currently not exposed

MOI.canget(instance::OSQPInstance, ::MOI.SolveTime) = hasresults(instance)
MOI.get(instance::OSQPInstance, ::MOI.SolveTime) = instance.results.info.run_time

MOI.canget(instance::OSQPInstance, ::MOI.TerminationStatus) = hasresults(instance)
function MOI.get(instance::OSQPInstance, ::MOI.TerminationStatus)
    osqpstatus = instance.results.info.status
    if osqpstatus == :Unsolved
        error("Problem is unsolved.") # TODO: good idea?
    elseif osqpstatus == :Interrupted
        MOI.Interrupted
    elseif osqpstatus == :Dual_infeasible
        MOI.InfeasibleOrUnbounded
    elseif osqpstatus == :Primal_infeasible
        MOI.InfeasibleNoResult
    elseif osqpstatus == :Max_iter_reached
        MOI.IterationLimit
    elseif osqpstatus == :Solved
        MOI.Success
    elseif osqpstatus == :Solved_inaccurate
        MOI.AlmostSuccess
    elseif osqpstatus == :Primal_infeasible_inaccurate
        MOI.InfeasibleNoResult
    end
end

MOI.canget(instance::OSQPInstance, ::MOI.PrimalStatus) = hasresults(instance)
function MOI.get(instance::OSQPInstance, ::MOI.PrimalStatus)
    osqpstatus = instance.results.info.status
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

MOI.canget(instance::OSQPInstance, ::MOI.DualStatus) = hasresults(instance)
function MOI.get(instance::OSQPInstance, ::MOI.DualStatus)
    osqpstatus = instance.results.info.status
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
MOI.candelete(instance::OSQPInstance, index::MOI.Index) = false
MOI.isvalid(instance::OSQPInstance, vi::VI) = vi.value ∈ 1 : get(instance, MOI.NumberOfVariables())
MOI.canaddvariable(instance::OSQPInstance) = false


## Variable attributes:
MOI.canget(instance::OSQPInstance, ::MOI.VariablePrimalStart, ::Type{VI}) = false # currently not exposed, but could be
MOI.canset(instance::OSQPInstance, ::MOI.VariablePrimalStart, ::Type{VI}) = false # TODO: need selective way of updating primal start

MOI.canget(instance::OSQPInstance, ::MOI.VariablePrimal, ::VI) = hasresults(instance) && instance.results.info.status ∈ OSQP.SOLUTION_PRESENT
MOI.get(instance::OSQPInstance, ::MOI.VariablePrimal, vi::VI) = instance.results.x[vi.value]
MOI.get(instance::OSQPInstance, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(instance, a, vi) # TODO: copied from SCS. Necessary?


## Constraints:
MOI.isvalid(instance::OSQPInstance, ci::CI) = ci.value ∈ 1 : get(instance, MOI.NumberOfConstraints())
MOI.canaddconstraint(instance::OSQPInstance, ::Type{F}, ::Type{S}) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet} = false

# TODO: can't modifyconstraint! with AbstractFunction because selective way to update A not exposed
MOI.canmodifyconstraint(instance::OSQPInstance, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, ::Type{MOI.ScalarAffineFunction}) = false
# MOI.modifyconstraint!(instance::OSQPInstance, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, func::MOI.ScalarAffineFunction)

# TODO: can't modifyconstraint! with AbstractSet because selective way to update lb and ub not exposed
MOI.canmodifyconstraint(instance::OSQPInstance, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, ::Type{MOI.Interval}) = false
# MOI.modifyconstraint!(instance::OSQPInstance, ci::CI{MOI.ScalarAffineFunction, MOI.Interval}, set::MOI.Interval)

# TODO: partial change with MultirowChange

MOI.supportsconstraint(instance::OSQPInstance, ::Type{MOI.ScalarAffineFunction}, ::Type{MOI.Interval}) = true


## Constraint attributes:
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintPrimalStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintDualStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canset(instance::OSQPInstance, ::MOI.ConstraintDualStart, ::Type{VI}) = false # TODO: need selective way of updating primal start
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintPrimal, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintDual, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintFunction, ::Type{<:CI}) = false # TODO
MOI.canget(instance::OSQPInstance, ::MOI.ConstraintSet, ::Type{<:CI}) = false # TODO


# Objective modification
MOI.canmodifyobjective(instance::OSQPInstance, ::Type{MOI.ScalarCoefficientChange}) = false # TODO: selective way of updating objective coefficients not exposed

end # module
