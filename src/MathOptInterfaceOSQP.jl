module MathOptInterfaceOSQP

export OSQPOptimizer, OSQPSettings

using Compat
using MathOptInterface
using MathOptInterfaceUtilities

const MOI = MathOptInterface
const MOIU = MathOptInterfaceUtilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SparseTriplets = Tuple{Vector{<:Integer}, Vector{<:Integer}, Vector{<:Any}}

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

    MOI.empty!(dest)
    idxmap = MOIU.IndexMap(dest, src)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        MOI.supportsconstraint(dest, F, S) || return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Unsupported $F-in-$S constraint", idxmap)
    end
    dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
    A, l, u = processconstraints(src, idxmap)
    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, dest.settings...)

    # TODO: clean up:
    # Process variable attributes
    m, n = size(A)
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

"""
Set up index map from `src` variables and constraints to `dest` variables and constraints.
"""
function MOIU.IndexMap(dest::OSQPOptimizer, src::MOI.ModelLike)
    idxmap = MOIU.IndexMap()
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = VI(i)
    end
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = CI{F, S}(i)
        end
    end
    idxmap
end

"""
Return objective sense, as well as matrix `P`, vector `q`, and scalar `c` such that objective function is 1/2 `x' P x + q' x + c`.
"""
function processobjective(src::MOI.ModelLike, idxmap)
    sense = MOI.get(src, MOI.ObjectiveSense())
    n = MOI.get(src, MOI.NumberOfVariables())

    if sense != MOI.FeasibilitySense
        if MOI.canget(src, MOI.ObjectiveFunction{SingleVariable}())
            fsingle = MOI.get(src, MOI.ObjectiveFunction{SingleVariable}())
            P = spzeros(n, n)
            q = zeros(n)
            q[idxmap[fsingle.variable.value]] = 1
            c = 0.
        elseif MOI.canget(src, MOI.ObjectiveFunction{Affine}())
            faffine = MOI.get(src, MOI.ObjectiveFunction{Affine}())
            P = spzeros(n, n)
            q = processlinearterm(faffine.variables, faffine.coefficients, idxmap)
            c = faffine.constant
        elseif MOI.canget(src, MOI.ObjectiveFunction{Quadratic}())
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{Quadratic}())
            I = [idxmap[var].value for var in fquadratic.quadratic_rowvariables]
            J = [idxmap[var].value for var in fquadratic.quadratic_colvariables]
            V = fquadratic.quadratic_coefficients
            symmetrize!((I, J, V))
            P = sparse(I, J, V, n, n)
            q = processlinearterm(fquadratic.affine_variables, fquadratic.affine_coefficients, idxmap)
            c = fquadratic.constant
        else
            error("No suitable objective function found")
        end
        sense == MOI.MaxSense && (scale!(P, -1); scale!(q, -1); c = -c)
    else
        P = spzeros(n, n)
        q = zeros(n)
        c = 0.
    end
    sense, P, q, c
end

function processlinearterm(variables::Vector{VI}, coefficients::Vector, idxmap)
    q = zeros(length(idxmap.varmap))
    ncoeffs = length(coefficients)
    @assert length(variables) == ncoeffs
    for i = 1 : ncoeffs
        q[idxmap[variables[i]].value] += coefficients[i]
    end
    q
end

function symmetrize!(triplets::SparseTriplets)
    I, J, V = triplets
    n = length(V)
    @assert length(I) == length(J) == n
    for i = 1 : n
        if I[i] != J[i]
            push!(I, J[i])
            push!(J, I[i])
            push!(V, V[i])
        end
    end
end

function processconstraints(src::MOI.ModelLike, idxmap)
    m = length(idxmap.conmap)
    l = Vector{Float64}(uninitialized, m)
    u = Vector{Float64}(uninitialized, m)
    bounds = (l, u)

    I = Int[]
    J = Int[]
    V = Float64[]
    triplets = (I, J, V)

    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!(triplets, bounds, src, idxmap, F, S)
    end
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)

    (A, l, u)
end

function processconstraints!(triplets::SparseTriplets, bounds::Tuple{<:Vector, <:Vector}, src::MOI.ModelLike, idxmap,
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    (l, u) = bounds
    for ci in cis_src
        s = MOI.Interval(MOI.get(src, MOI.ConstraintSet(), ci))
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        row = idxmap[ci].value
        l[row] = s.lower - constant(f)
        u[row] = s.upper - constant(f)
        processconstraints!(triplets, f, row, idxmap)
    end
end

function processconstraints!(triplets::SparseTriplets, f::MOI.SingleVariable, row::Int, idxmap)
    (I, J, V) = triplets
    col = idxmap[f.variable].value
    push!(I, row)
    push!(J, col)
    push!(V, 1)
    nothing
end

function processconstraints!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap)
    (I, J, V) = triplets
    ncoeff = length(f.coefficients)
    @assert length(f.variables) == ncoeff
    for i = 1 : ncoeff
        var = f.variables[i]
        coeff = f.coefficients[i]
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, coeff)
    end
end

## Standard optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = optimizer.sense

MOI.canget(optimizer::OSQPOptimizer, ::MOI.NumberOfVariables) = !optimizer.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
MOI.get(optimizer::OSQPOptimizer, ::MOI.NumberOfVariables) = OSQP.dimensions(optimizer.model)[1]

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = MOI.canget(optimizer, MOI.NumberOfVariables())
MOI.get(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = [VI(i) for i = 1 : get(optimizer, MOI.NumberOfVariables())] # TODO: support for UnitRange would be nice

# FIXME or remove:
# MOI.canget(optimizer::OSQPOptimizer, ::MOI.NumberOfConstraints) = !optimizer.isempty # https://github.com/oxfordcontrol/OSQP.jl/issues/10
# MOI.get(optimizer::OSQPOptimizer, ::MOI.NumberOfConstraints) = OSQP.dimensions(optimizer.model)[2]

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
# MOI.isvalid(optimizer::OSQPOptimizer, ci::CI) = ci.value ∈ 1 : get(optimizer, MOI.NumberOfConstraints()) # FIXME
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
