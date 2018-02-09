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

    # Process instance attributes
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

function processconstraints!(A::SparseMatrixCSC, l::Vector, u::Vector, src::MOI.AbstractInstance, idxmap, F::Type{<:MOI.ScalarAffineFunction}, S::Type{<:MOI.Interval})
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

MOI.canget(instance::OSQPInstance, ::MOI.VariablePrimal, ::Type{VI}) = hasresults(instance) && instance.results.info.status ∈ OSQP.SOLUTION_PRESENT
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

MOI.supportsconstraint(instance::OSQPInstance, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true


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
