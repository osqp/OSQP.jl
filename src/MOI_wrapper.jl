module MathOptInterfaceOSQP

include("modcaches.jl")
using .ModificationCaches

import LinearAlgebra: rmul!
import MathOptInterface
import SparseArrays
import OSQP

export Optimizer, OSQPSettings, OSQPModel

const MOI = MathOptInterface

const SparseTriplets = Tuple{Vector{Int},Vector{Int},Vector{<:Any}}

const IntervalConvertible = Union{
    MOI.Interval{Float64},
    MOI.LessThan{Float64},
    MOI.GreaterThan{Float64},
    MOI.EqualTo{Float64},
}

const SupportedVectorSets = Union{MOI.Zeros,MOI.Nonnegatives,MOI.Nonpositives}

lower(::MOI.Zeros, i::Int) = 0.0
lower(::MOI.Nonnegatives, i::Int) = 0.0
lower(::MOI.Nonpositives, i::Int) = -Inf
upper(::MOI.Zeros, i::Int) = 0.0
upper(::MOI.Nonnegatives, i::Int) = Inf
upper(::MOI.Nonpositives, i::Int) = 0.0

# TODO: just use ∈ on 0.7 (allocates on 0.6):
function _contains(haystack, needle)
    for x in haystack
        x == needle && return true
    end
    return false
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::OSQP.Model
    hasresults::Bool
    results::OSQP.Results
    is_empty::Bool
    silent::Bool
    settings::Dict{Symbol,Any} # need to store these, because they should be preserved if empty! is called
    sense::MOI.OptimizationSense
    objconstant::Float64
    constrconstant::Vector{Float64}
    modcache::ProblemModificationCache{Float64}
    warmstartcache::WarmStartCache{Float64}
    rowranges::Dict{Int,UnitRange{Int}}

    function Optimizer(; kwargs...)
        inner = OSQP.Model()
        hasresults = false
        results = OSQP.Results()
        is_empty = true
        sense = MOI.MIN_SENSE
        objconstant = 0.0
        constrconstant = Float64[]
        modcache = ProblemModificationCache{Float64}()
        warmstartcache = WarmStartCache{Float64}()
        rowranges = Dict{Int,UnitRange{Int}}()
        optimizer = new(
            inner,
            hasresults,
            results,
            is_empty,
            false,
            Dict{Symbol,Any}(:verbose => true),
            sense,
            objconstant,
            constrconstant,
            modcache,
            warmstartcache,
            rowranges,
        )
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawOptimizerAttribute(key), value)
        end
        return optimizer
    end
end

MOI.get(::Optimizer, ::MOI.SolverName) = "OSQP"

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
    if !MOI.is_empty(optimizer)
        if optimizer.silent
            OSQP.update_settings!(optimizer.inner; :verbose => false)
        else
            OSQP.update_settings!(
                optimizer.inner;
                :verbose => optimizer.settings[:verbose],
            )
        end
    end
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, limit::Real)
    MOI.set(model, OSQPSettings.TimeLimit(), limit)
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(model.settings, :time_limit)
    if !MOI.is_empty(model)
        OSQP.update_settings!(model.inner, time_limit = 0.0)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.settings, :time_limit, nothing)
end

hasresults(optimizer::Optimizer) = optimizer.hasresults

function MOI.empty!(optimizer::Optimizer)
    optimizer.inner = OSQP.Model()
    optimizer.hasresults = false
    optimizer.results = OSQP.Results()
    optimizer.sense = MOI.MIN_SENSE # model parameter, so needs to be reset
    optimizer.objconstant = 0.0
    optimizer.constrconstant = Float64[]
    optimizer.modcache = ProblemModificationCache{Float64}()
    optimizer.warmstartcache = WarmStartCache{Float64}()
    empty!(optimizer.rowranges)
    optimizer
    return
end

MOI.is_empty(optimizer::Optimizer) = optimizer.inner.isempty

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    idxmap = MOI.Utilities.IndexMap(dest, src)
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
    A, l, u, dest.constrconstant =
        processconstraints(src, idxmap, dest.rowranges)
    settings = copy(dest.settings)
    if dest.silent
        settings[:verbose] = false
    end
    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, settings...)
    dest.modcache = ProblemModificationCache(P, q, A, l, u)
    dest.warmstartcache = WarmStartCache{Float64}(size(A, 2), size(A, 1))
    processprimalstart!(dest.warmstartcache.x, src, idxmap)
    processdualstart!(dest.warmstartcache.y, src, idxmap, dest.rowranges)
    return idxmap
end

"""
Set up index map from `src` variables and constraints to `dest` variables and
constraints.
"""
function MOI.Utilities.IndexMap(dest::Optimizer, src::MOI.ModelLike)
    idxmap = MOI.Utilities.IndexMap()
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.supports_constraint(dest, F, S) ||
            throw(MOI.UnsupportedConstraint{F,S}())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = MOI.ConstraintIndex{F,S}(i)
        end
    end
    return idxmap
end

function assign_constraint_row_ranges!(
    rowranges::Dict{Int,UnitRange{Int}},
    idxmap::MOI.Utilities.IndexMap,
    src::MOI.ModelLike,
)
    startrow = 1
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        for ci_src in cis_src
            set = MOI.get(src, MOI.ConstraintSet(), ci_src)
            ci_dest = idxmap[ci_src]
            endrow = startrow + MOI.dimension(set) - 1
            rowranges[ci_dest.value] = startrow:endrow
            startrow = endrow + 1
        end
    end
    return
end

function constraint_rows(
    rowranges::Dict{Int,UnitRange{Int}},
    ci::MOI.ConstraintIndex{<:Any,<:MOI.AbstractScalarSet},
)
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    return first(rowrange)
end

function constraint_rows(
    rowranges::Dict{Int,UnitRange{Int}},
    ci::MOI.ConstraintIndex{<:Any,<:MOI.AbstractVectorSet},
)
    return rowranges[ci.value]
end

function constraint_rows(optimizer::Optimizer, ci::MOI.ConstraintIndex)
    return constraint_rows(optimizer.rowranges, ci)
end

"""
Return objective sense, as well as matrix `P`, vector `q`, and scalar `c` such
that objective function is `1/2 x' P x + q' x + c`.
"""
function processobjective(src::MOI.ModelLike, idxmap)
    sense = MOI.get(src, MOI.ObjectiveSense())
    n = MOI.get(src, MOI.NumberOfVariables())
    q = zeros(n)
    if sense != MOI.FEASIBILITY_SENSE
        function_type = MOI.get(src, MOI.ObjectiveFunctionType())
        if function_type == MOI.ScalarAffineFunction{Float64}
            faffine = MOI.get(
                src,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            )
            P = SparseArrays.spzeros(n, n)
            processlinearterms!(q, faffine.terms, idxmap)
            c = faffine.constant
        elseif function_type == MOI.ScalarQuadraticFunction{Float64}
            fquadratic = MOI.get(
                src,
                MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(),
            )
            I = [
                Int(idxmap[term.variable_1].value) for
                term in fquadratic.quadratic_terms
            ]
            J = [
                Int(idxmap[term.variable_2].value) for
                term in fquadratic.quadratic_terms
            ]
            V = [term.coefficient for term in fquadratic.quadratic_terms]
            upper_triangularize!(I, J, V)
            P = SparseArrays.sparse(I, J, V, n, n)
            processlinearterms!(q, fquadratic.affine_terms, idxmap)
            c = fquadratic.constant
        else
            throw(
                MOI.UnsupportedAttribute(
                    MOI.ObjectiveFunction{function_type}(),
                ),
            )
        end
        sense == MOI.MAX_SENSE && (rmul!(P, -1); rmul!(q, -1); c = -c)
    else
        P = SparseArrays.spzeros(n, n)
        q = zeros(n)
        c = 0.0
    end
    return sense, P, q, c
end

function processlinearterms!(
    q,
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmapfun::Function = identity,
)
    # This is currently needed to avoid depwarns. TODO: make this nice again:
    if q isa VectorModificationCache
        q[:] = 0
    else
        q .= 0
    end
    for term in terms
        var = term.variable
        coeff = term.coefficient
        q[idxmapfun(var).value] += coeff
    end
    return
end

function processlinearterms!(
    q,
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmap::MOI.Utilities.IndexMap,
)
    processlinearterms!(q, terms, var -> idxmap[var])
    return
end

function upper_triangularize!(I::Vector{Int}, J::Vector{Int}, V::Vector)
    n = length(V)
    (length(I) == length(J) == n) || error()
    for i in 1:n
        if I[i] > J[i]
            I[i], J[i] = J[i], I[i]
        end
    end
    return
end

function processconstraints(
    src::MOI.ModelLike,
    idxmap,
    rowranges::Dict{Int,UnitRange{Int}},
)
    m = mapreduce(length, +, values(rowranges), init = 0)
    l = Vector{Float64}(undef, m)
    u = Vector{Float64}(undef, m)
    constant = Vector{Float64}(undef, m)
    bounds = (l, u)
    I = Int[]
    J = Int[]
    V = Float64[]
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        processconstraints!(
            (I, J, V),
            bounds,
            constant,
            src,
            idxmap,
            rowranges,
            F,
            S,
        )
    end
    l .-= constant
    u .-= constant
    n = MOI.get(src, MOI.NumberOfVariables())
    A = SparseArrays.sparse(I, J, V, m, n)
    return A, l, u, constant
end

function processconstraints!(
    triplets::SparseTriplets,
    bounds::Tuple{<:Vector,<:Vector},
    constant::Vector{Float64},
    src::MOI.ModelLike,
    idxmap,
    rowranges::Dict{Int,UnitRange{Int}},
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_rows(rowranges, idxmap[ci])
        processconstant!(constant, rows, f)
        processlinearpart!(triplets, f, rows, idxmap)
        processconstraintset!(bounds, rows, s)
    end
    return
end

function processconstant!(
    c::Vector{Float64},
    row::Int,
    f::MOI.ScalarAffineFunction{Float64},
)
    c[row] = MOI.constant(f, Float64)
    return
end

function processconstant!(
    c::Vector{Float64},
    rows::UnitRange{Int},
    f::MOI.VectorAffineFunction{Float64},
)
    for (i, row) in enumerate(rows)
        c[row] = f.constants[i]
    end
    return
end

function processlinearpart!(
    triplets::SparseTriplets,
    f::MOI.ScalarAffineFunction,
    row::Int,
    idxmap,
)
    (I, J, V) = triplets
    for term in f.terms
        var = term.variable
        coeff = term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, coeff)
    end
    return
end

function processlinearpart!(
    triplets::SparseTriplets,
    f::MOI.VectorAffineFunction,
    rows::UnitRange{Int},
    idxmap,
)
    (I, J, V) = triplets
    for term in f.terms
        row = rows[term.output_index]
        var = term.scalar_term.variable
        coeff = term.scalar_term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, coeff)
    end
    return
end

function processconstraintset!(
    bounds::Tuple{<:Vector,<:Vector},
    row::Int,
    s::IntervalConvertible,
)
    processconstraintset!(bounds, row, MOI.Interval(s))
    return
end

function processconstraintset!(
    bounds::Tuple{<:Vector,<:Vector},
    row::Int,
    interval::MOI.Interval{Float64},
)
    l, u = bounds
    l[row] = interval.lower
    u[row] = interval.upper
    return
end

function processconstraintset!(
    bounds::Tuple{<:Vector,<:Vector},
    rows::UnitRange{Int},
    s::S,
) where {S<:SupportedVectorSets}
    l, u = bounds
    for (i, row) in enumerate(rows)
        l[row] = lower(s, i)
        u[row] = upper(s, i)
    end
    return
end

function processprimalstart!(x, src::MOI.ModelLike, idxmap)
    has_primal_start = false
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr isa MOI.VariablePrimalStart
            has_primal_start = true
        end
    end
    if has_primal_start
        vis_src = MOI.get(src, MOI.ListOfVariableIndices())
        for vi in vis_src
            value = MOI.get(src, MOI.VariablePrimalStart(), vi)
            if value !== nothing
                x[idxmap[vi].value] = value
            end
        end
    end
    return
end

function processdualstart!(
    y,
    src::MOI.ModelLike,
    idxmap,
    rowranges::Dict{Int,UnitRange{Int}},
)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        has_dual_start = false
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
            if attr isa MOI.ConstraintDualStart
                has_dual_start = true
            end
        end
        if has_dual_start
            cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
            for ci in cis_src
                rows = constraint_rows(rowranges, idxmap[ci])
                dual = MOI.get(src, MOI.ConstraintDualStart(), ci)
                if dual !== nothing
                    for (i, row) in enumerate(rows)
                        y[row] = -dual[i] # opposite dual convention
                    end
                end
            end
        end
    end
    return
end

## Standard optimizer attributes:

MOI.get(optimizer::Optimizer, ::MOI.ObjectiveSense) = optimizer.sense

function MOI.get(optimizer::Optimizer, ::MOI.NumberOfVariables)
    return OSQP.dimensions(optimizer.inner)[1]
end

function MOI.get(optimizer::Optimizer, ::MOI.ListOfVariableIndices)
    return [
        MOI.VariableIndex(i) for
        i in 1:MOI.get(optimizer, MOI.NumberOfVariables())
    ] # TODO: support for UnitRange would be nice
end

## Solver-specific optimizer attributes:
module OSQPSettings

using MathOptInterface
using OSQP

export OSQPAttribute, isupdatable

abstract type OSQPAttribute <: MathOptInterface.AbstractOptimizerAttribute end

# TODO: just use ∈ on 0.7 (allocates on 0.6):
function _contains(haystack, needle)
    for x in haystack
        x == needle && return true
    end
    return false
end

for setting in fieldnames(OSQP.Settings)
    Attribute =
        Symbol(mapreduce(uppercasefirst, *, split(String(setting), '_'))) # to camelcase
    @eval begin
        export $Attribute
        struct $Attribute <: OSQPAttribute end
        Base.Symbol(::$Attribute) = $(QuoteNode(setting))
        function isupdatable(::$Attribute)
            return $(_contains(OSQP.UPDATABLE_SETTINGS, setting))
        end
    end
end
end # module

using .OSQPSettings

_symbol(param::MOI.RawOptimizerAttribute) = Symbol(param.name)

_symbol(a::OSQPAttribute) = Symbol(a)

function OSQPSettings.isupdatable(param::MOI.RawOptimizerAttribute)
    return _contains(OSQP.UPDATABLE_SETTINGS, _symbol(param))
end

function MOI.set(
    optimizer::Optimizer,
    a::Union{OSQPAttribute,MOI.RawOptimizerAttribute},
    value,
)
    (isupdatable(a) || MOI.is_empty(optimizer)) ||
        throw(MOI.SetAttributeNotAllowed(a))
    setting = _symbol(a)
    optimizer.settings[setting] = value
    if !MOI.is_empty(optimizer)
        OSQP.update_settings!(optimizer.inner; setting => value)
    end
end

function MOI.get(
    optimizer::Optimizer,
    a::Union{OSQPAttribute,MOI.RawOptimizerAttribute},
)
    return optimizer.settings[_symbol(a)]
end

## Optimizer methods:
function MOI.optimize!(optimizer::Optimizer)
    processupdates!(optimizer.inner, optimizer.modcache)
    processupdates!(optimizer.inner, optimizer.warmstartcache)
    OSQP.solve!(optimizer.inner, optimizer.results)
    optimizer.hasresults = true
    # Copy previous solution into warm start cache without setting the dirty
    # bit:
    copyto!(optimizer.warmstartcache.x.data, optimizer.results.x)
    copyto!(optimizer.warmstartcache.y.data, optimizer.results.y)
    return
end

## Optimizer attributes:

MOI.get(optimizer::Optimizer, ::MOI.RawSolver) = optimizer.inner

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = optimizer.hasresults ? 1 : 0

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    optimizer::Optimizer,
    a::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    obj::MOI.ScalarAffineFunction{Float64},
)
    if MOI.is_empty(optimizer)
        throw(MOI.SetAttributeNotAllowed(a))
    end
    optimizer.modcache.P[:] = 0
    processlinearterms!(optimizer.modcache.q, obj.terms)
    optimizer.objconstant = MOI.constant(obj)
    return
end

function MOI.set(
    optimizer::Optimizer,
    a::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    obj::MOI.ScalarQuadraticFunction{Float64},
)
    MOI.is_empty(optimizer) && throw(MOI.SetAttributeNotAllowed(a))
    cache = optimizer.modcache
    cache.P[:] = 0
    for term in obj.quadratic_terms
        row = term.variable_1.value
        col = term.variable_2.value
        coeff = term.coefficient
        row > col && ((row, col) = (col, row)) # upper triangle only
        if !(CartesianIndex(row, col) in cache.P.cartesian_indices_set)
            throw(
                MOI.SetAttributeNotAllowed(
                    a,
                    "This nonzero entry was not in the sparsity pattern of " *
                    "the objective function provided at `MOI.copy_to` and " *
                    "OSQP does not support changing the sparsity pattern.",
                ),
            )
        end
        cache.P[row, col] += coeff
    end
    processlinearterms!(optimizer.modcache.q, obj.affine_terms)
    optimizer.objconstant = MOI.constant(obj)
    return
end

function MOI.get(optimizer::Optimizer, a::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, a)
    rawobj = optimizer.results.info.obj_val + optimizer.objconstant
    return ifelse(optimizer.sense == MOI.MAX_SENSE, -rawobj, rawobj)
end

function check_has_results(optimizer::Optimizer)
    if !hasresults(optimizer)
        error("Problem is unsolved.")
    end
    return
end

# Since these aren't explicitly returned by OSQP, I feel like it would be better
# to have a fallback method compute these:
function MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec)
    check_has_results(optimizer)
    return optimizer.results.info.run_time
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return string(optimizer.results.info.status)
end

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    hasresults(optimizer) || return MOI.OPTIMIZE_NOT_CALLED
    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        return MOI.OPTIMIZE_NOT_CALLED
    elseif osqpstatus == :Interrupted
        return MOI.INTERRUPTED
    elseif osqpstatus == :Dual_infeasible
        return MOI.DUAL_INFEASIBLE
    elseif osqpstatus == :Primal_infeasible
        return MOI.INFEASIBLE
    elseif osqpstatus == :Max_iter_reached
        return MOI.ITERATION_LIMIT
    elseif osqpstatus == :Solved
        return MOI.OPTIMAL
    elseif osqpstatus == :Solved_inaccurate
        return MOI.ALMOST_OPTIMAL
    elseif osqpstatus == :Primal_infeasible_inaccurate
        return MOI.ALMOST_INFEASIBLE
    else
        @assert osqpstatus == :Non_convex
        return MOI.INVALID_MODEL
    end
end

function MOI.get(optimizer::Optimizer, a::MOI.PrimalStatus)
    if a.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        return MOI.NO_SOLUTION
    elseif osqpstatus == :Primal_infeasible
        # FIXME is it `NO_SOLUTION` (e.g. `NaN`s) or `INFEASIBLE_POINT` (e.g.
        # current primal solution that we know is infeasible with the dual
        # certificate)
        return MOI.NO_SOLUTION
    elseif osqpstatus == :Solved
        return MOI.FEASIBLE_POINT
    elseif osqpstatus == :Primal_infeasible_inaccurate
        return MOI.UNKNOWN_RESULT_STATUS
    elseif osqpstatus == :Dual_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex
        # (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        return MOI.NO_SOLUTION
    end
end

function MOI.get(optimizer::Optimizer, a::MOI.DualStatus)
    if a.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    osqpstatus = optimizer.results.info.status
    if osqpstatus == :Unsolved
        return MOI.NO_SOLUTION
    elseif osqpstatus == :Dual_infeasible
        # FIXME is it `NO_SOLUTION` (e.g. `NaN`s) or `INFEASIBLE_POINT` (e.g.
        # current dual solution that we know is infeasible with the dual
        # certificate)
        return MOI.NO_SOLUTION
    elseif osqpstatus == :Primal_infeasible
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif osqpstatus == :Primal_infeasible_inaccurate
        return MOI.NEARLY_INFEASIBILITY_CERTIFICATE
    elseif osqpstatus == :Solved
        return MOI.FEASIBLE_POINT
    else
        # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex
        # (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        return MOI.NO_SOLUTION
    end
end

## Variables:

function MOI.is_valid(optimizer::Optimizer, vi::MOI.VariableIndex)
    return vi.value ∈ 1:MOI.get(optimizer, MOI.NumberOfVariables())
end

## Variable attributes:

function MOI.get(
    optimizer::Optimizer,
    a::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, a)
    x = ifelse(
        _contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status),
        optimizer.results.x,
        optimizer.results.dual_inf_cert,
    )
    return x[vi.value]
end

function MOI.set(
    optimizer::Optimizer,
    a::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value,
)
    if MOI.is_empty(optimizer)
        throw(MOI.SetAttributeNotAllowed(a))
    end
    return optimizer.warmstartcache.x[vi.value] = value
end

## Constraints:

function MOI.is_valid(optimizer::Optimizer, ci::MOI.ConstraintIndex)
    if MOI.is_empty(optimizer)
        return false
    end
    return ci.value ∈ keys(optimizer.rowranges)
end

function MOI.set(
    optimizer::Optimizer,
    a::MOI.ConstraintDualStart,
    ci::MOI.ConstraintIndex,
    value,
)
    MOI.is_empty(optimizer) && throw(MOI.SetAttributeNotAllowed(a))
    rows = constraint_rows(optimizer, ci)
    for (i, row) in enumerate(rows)
        optimizer.warmstartcache.y[row] = -value[i] # opposite dual convention
    end
    return
end

# function modification:

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{
        MOI.ScalarAffineFunction{Float64},
        <:IntervalConvertible,
    },
    f::MOI.ScalarAffineFunction{Float64},
)
    if !MOI.is_valid(optimizer, ci)
        throw(MOI.InvalidIndex(ci))
    end
    row = constraint_rows(optimizer, ci)
    optimizer.modcache.A[row, :] = 0
    for term in f.terms
        col = term.variable.value
        coeff = term.coefficient
        optimizer.modcache.A[row, col] += coeff
    end
    Δconstant = optimizer.constrconstant[row] - f.constant
    optimizer.constrconstant[row] = f.constant
    optimizer.modcache.l[row] += Δconstant
    optimizer.modcache.u[row] += Δconstant
    return
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{
        MOI.VectorAffineFunction{Float64},
        <:SupportedVectorSets,
    },
    f::MOI.VectorAffineFunction{Float64},
)
    if !MOI.is_valid(optimizer, ci)
        throw(MOI.InvalidIndex(ci))
    end
    rows = constraint_rows(optimizer, ci)
    for row in rows
        optimizer.modcache.A[row, :] = 0
    end
    for term in f.terms
        row = rows[term.output_index]
        col = term.scalar_term.variable.value
        coeff = term.scalar_term.coefficient
        optimizer.modcache.A[row, col] += coeff
    end
    for (i, row) in enumerate(rows)
        Δconstant = optimizer.constrconstant[row] - f.constants[i]
        optimizer.constrconstant[row] = f.constants[i]
        optimizer.modcache.l[row] += Δconstant
        optimizer.modcache.u[row] += Δconstant
    end
    return
end

# set modification:

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
    s::S,
) where {S<:IntervalConvertible}
    if !MOI.is_valid(optimizer, ci)
        throw(MOI.InvalidIndex(ci))
    end
    interval = S <: MOI.Interval ? s : MOI.Interval(s)
    row = constraint_rows(optimizer, ci)
    constant = optimizer.constrconstant[row]
    optimizer.modcache.l[row] = interval.lower - constant
    optimizer.modcache.u[row] = interval.upper - constant
    return
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
    s::S,
) where {S<:SupportedVectorSets}
    if !MOI.is_valid(optimizer, ci)
        throw(MOI.InvalidIndex(ci))
    end
    rows = constraint_rows(optimizer, ci)
    for (i, row) in enumerate(rows)
        constant = optimizer.constrconstant[row]
        optimizer.modcache.l[row] = lower(s, i) - constant
        optimizer.modcache.u[row] = upper(s, i) - constant
    end
    return
end

# partial function modification:

function MOI.modify(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex{
        MOI.ScalarAffineFunction{Float64},
        <:IntervalConvertible,
    },
    change::MOI.ScalarCoefficientChange,
)
    MOI.is_valid(optimizer, ci) || throw(MOI.InvalidIndex(ci))
    row = constraint_rows(optimizer, ci)
    optimizer.modcache.A[row, change.variable.value] = change.new_coefficient
    return
end

# TODO: MultirowChange?

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:IntervalConvertible},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{Float64}},
    ::Type{<:SupportedVectorSets},
)
    return true
end

## Constraint attributes:

function MOI.get(
    optimizer::Optimizer,
    a::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex,
)
    MOI.check_result_index_bounds(optimizer, a)
    y = ifelse(
        _contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status),
        optimizer.results.y,
        optimizer.results.prim_inf_cert,
    )
    rows = constraint_rows(optimizer, ci)
    return -y[rows]
end

# Objective modification

function MOI.modify(
    optimizer::Optimizer,
    ::MOI.ObjectiveFunction,
    change::MOI.ScalarConstantChange,
)
    if MOI.is_empty(optimizer)
        throw(MOI.ModifyObjectiveNotAllowed(change))
    end
    constant = change.new_constant
    if optimizer.sense == MOI.MAX_SENSE
        constant = -constant
    end
    return optimizer.objconstant = constant
end

function MOI.modify(
    optimizer::Optimizer,
    ::MOI.ObjectiveFunction,
    change::MOI.ScalarCoefficientChange,
)
    if MOI.is_empty(optimizer)
        throw(MOI.ModifyObjectiveNotAllowed(change))
    end
    coef = change.new_coefficient
    if optimizer.sense == MOI.MAX_SENSE
        coef = -coef
    end
    return optimizer.modcache.q[change.variable.value] = coef
end

# There is currently no ScalarQuadraticCoefficientChange.

MOI.Utilities.@model(
    OSQPModel,
    (),
    (MOI.Interval, MOI.LessThan, MOI.GreaterThan, MOI.EqualTo),
    (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives),
    (),
    (),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (),
    (MOI.VectorAffineFunction,)
)

end # module
