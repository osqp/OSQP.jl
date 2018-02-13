module MathOptInterfaceOSQP

export OSQPOptimizer, OSQPSettings

using Compat
using MathOptInterface
using MathOptInterfaceUtilities

const MOI = MathOptInterface
const MOIU = MathOptInterfaceUtilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SingleVariable = MOI.SingleVariable
const Affine = MOI.ScalarAffineFunction{Float64}
const Quadratic = MOI.ScalarQuadraticFunction{Float64}
const AffineConvertible = Union{Affine, SingleVariable}

const Interval = MOI.Interval{Float64}
const LessThan = MOI.LessThan{Float64}
const GreaterThan = MOI.GreaterThan{Float64}
const EqualTo = MOI.EqualTo{Float64}
const IntervalConvertible = Union{Interval, LessThan, GreaterThan, EqualTo}

import OSQP
import MathOptInterfaceUtilities: IndexMap

# TODO: consider moving to MOI:
constant(f::SingleVariable) = 0
constant(f::Affine) = f.constant
constant(f::Quadratic) = f.constant

mutable struct OSQPOptimizer <: MOI.AbstractOptimizer # TODO: name?
    inner::OSQP.Model
    results::Union{OSQP.Results, Nothing}
    isempty::Bool
    settings::Dict{Symbol, Any} # need to store these, because they should be preserved if empty! is called
    sense::MOI.OptimizationSense
    objconstant::Float64

    OSQPOptimizer() = new(OSQP.Model(), nothing, true, Dict{Symbol, Any}(), MOI.MinSense, 0.)
end

hasresults(optimizer::OSQPOptimizer) = optimizer.results != nothing

function MOI.empty!(optimizer::OSQPOptimizer)
    optimizer.inner = OSQP.Model()
    optimizer.results = nothing
    optimizer.isempty = true
    optimizer.sense = MOI.MinSense # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer
end

MOI.isempty(optimizer::OSQPOptimizer) = optimizer.isempty

function MOI.copy!(dest::OSQPOptimizer, src::MOI.ModelLike)
    # TODO: cangets (awaiting https://github.com/JuliaOpt/MathOptInterface.jl/pull/210)

    # Empty
    MOI.empty!(dest)

    # Compute problem size, set up index map, and check constraint types
    n = MOI.get(src, MOI.NumberOfVariables())
    idxmap = IndexMap()
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = VI(i)
    end
    m = 0
    let i = 0
        for (F, S) in MOI.get(src, MOI.ListOfConstraints())
            MOI.supportsconstraint(dest, F, S) || return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Unsupported $F-in-$S constraint", idxmap)
            m += MOI.get(src, MOI.NumberOfConstraints{F, S}())
            cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            for ci in cis_src
                i += 1
                idxmap[ci] = CI{F, S}(i)
            end
        end
    end

    # Allocate storage for problem data
    P = spzeros(n, n)
    q = zeros(n)
    A = spzeros(m, n)
    l = zeros(m)
    u = zeros(m)

    # Process objective function
    # prefer the simplest form
    sense = dest.sense = MOI.get(src, MOI.ObjectiveSense())
    if sense != MOI.FeasibilitySense
        if MOI.canget(src, MOI.ObjectiveFunction{SingleVariable}())
            objfun = MOI.get(src, MOI.ObjectiveFunction{SingleVariable}())
            processobjective!(P, q, objfun, idxmap)
        elseif MOI.canget(src, MOI.ObjectiveFunction{Affine}())
            objfun = MOI.get(src, MOI.ObjectiveFunction{Affine}())
            processobjective!(P, q, objfun, idxmap)
            dest.objconstant = constant(objfun)
        elseif MOI.canget(src, MOI.ObjectiveFunction{Quadratic}())
            objfun = MOI.get(src, MOI.ObjectiveFunction{Quadratic}())
            processobjective!(P, q, objfun, idxmap)
            dest.objconstant = constant(objfun)
        else
            return MOI.CopyResult(MOI.CopyOtherError, "No suitable objective function found", idxmap)
        end
        sense == MOI.MaxSense && (scale!(P, -1); scale!(q, -1); dest.objconstant = -dest.objconstant)
    end

    # Process constraints
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!(A, l, u, src, idxmap, F, S)
    end

    # Load data into OSQP Model
    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, dest.settings...)

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
        for (F, S) in MOI.get(src, MOI.ListOfConstraints())
            cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            for ci in cis_src
                y[idxmap[ci].value] = get(src, MOI.ConstraintDualStart(), ci)
            end
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
        q[idxmap[variables[i]].value] += coefficients[i]
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

function processconstraintfun!(A::AbstractMatrix, row::Int, idxmap, f::SingleVariable)
    col = idxmap[f.variable].value
    A[row, col] = 1
    nothing
end

function processconstraintfun!(A::AbstractMatrix, row::Int, idxmap, f::Affine)
    n = length(f.coefficients)
    @assert length(f.variables) == n
    for i = 1 : n
        var = f.variables[i]
        coeff = f.coefficients[i]
        col = idxmap[var].value
        A[row, col] = coeff
    end
    nothing
end

function processconstraints!(A::SparseMatrixCSC, l::Vector, u::Vector, src::MOI.ModelLike, idxmap, F::Type, S::Type)
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    for ci in cis_src
        s = MOI.Interval(MOI.get(src, MOI.ConstraintSet(), ci))
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        row = idxmap[ci].value
        l[row] = s.lower - constant(f)
        u[row] = s.upper - constant(f)
        processconstraintfun!(A, row, idxmap, f)
    end
end

## Standard optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = optimizer.sense

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


## Solver-specific optimizer attributes:
module OSQPSettings

export OSQPAttribute, isupdatable

using MathOptInterface, OSQP
const MOI = MathOptInterface

abstract type OSQPAttribute <: MOI.AbstractOptimizerAttribute end

for setting in fieldnames(OSQP.Settings)
    Attribute = Symbol(mapreduce(ucfirst, *, split(String(setting), '_'))) # to camelcase
    @eval begin
        export $Attribute
        struct $Attribute <: OSQPAttribute end
        Base.Symbol(::$Attribute) = $(QuoteNode(setting))
        isupdatable(::$Attribute) = $(setting ∈ OSQP.UPDATABLE_SETTINGS)
    end
end
end # module

using .OSQPSettings

MOI.canset(optimizer::OSQPOptimizer, a::OSQPAttribute) = MOI.isempty(optimizer) || isupdatable(a)
function MOI.set!(optimizer::OSQPOptimizer, a::OSQPAttribute, value)
    MOI.canset(optimizer, a) || error()
    setting = Symbol(a)
    optimizer.settings[setting] = value
    if !MOI.isempty(optimizer)
        OSQP.update_settings!(optimizer.inner; setting = value)
    end
end

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
function MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveValue)
    rawobj = optimizer.results.info.obj_val + optimizer.objconstant
    ifelse(optimizer.sense == MOI.MaxSense, -rawobj, rawobj)
end

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
    elseif osqpstatus == :Dual_infeasible
        MOI.InfeasibilityCertificate
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
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
    elseif osqpstatus == :Primal_infeasible
        MOI.InfeasibilityCertificate
    elseif osqpstatus == :Primal_infeasible_inaccurate
        MOI.AlmostInfeasibilityCertificate
    elseif osqpstatus == :Solved
        MOI.FeasiblePoint
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
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

MOI.canget(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, ::Type{VI}) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, vi::VI)
    if optimizer.results.info.status in OSQP.SOLUTION_PRESENT
        return optimizer.results.x[vi.value]
    else
        if optimizer.results.dual_inf_cert != nothing
            return optimizer.results.dual_inf_cert[vi.value]
        else
            error("Variable primal not available")
        end
    end
end


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

MOI.supportsconstraint(optimizer::OSQPOptimizer, ::Type{<:AffineConvertible}, ::Type{<:IntervalConvertible}) = true


## Constraint attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintPrimalStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canset(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{VI}) = false # TODO: need selective way of updating primal start
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintPrimal, ::Type{<:CI}) = false # currently not exposed, but could be
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintFunction, ::Type{<:CI}) = false # TODO
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintSet, ::Type{<:CI}) = false # TODO

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDual, ::Type{<:CI}) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, ::MOI.ConstraintDual, ci::CI)
    # MOI uses opposite dual convention
    if optimizer.results.info.status in OSQP.SOLUTION_PRESENT
        return -optimizer.results.y[ci.value]
    else
        if optimizer.results.prim_inf_cert != nothing
            return -optimizer.results.prim_inf_cert[ci.value]
        else
            error("Constraint dual not available.")
        end
    end
end


# Objective modification
MOI.canmodifyobjective(optimizer::OSQPOptimizer, ::Type{MOI.ScalarCoefficientChange}) = false # TODO: selective way of updating objective coefficients not exposed

end # module
