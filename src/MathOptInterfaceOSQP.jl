module MathOptInterfaceOSQP

include("modcaches.jl")
using .ModificationCaches

using Compat
using Compat.SparseArrays
using MathOptInterface
using MathOptInterface.Utilities

export OSQPOptimizer, OSQPSettings, OSQPModel

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SparseTriplets = Tuple{Vector{Int}, Vector{Int}, Vector{<:Any}}

const SingleVariable = MOI.SingleVariable
const Affine = MOI.ScalarAffineFunction{Float64}
const Quadratic = MOI.ScalarQuadraticFunction{Float64}
const AffineConvertible = Union{Affine, SingleVariable}
const VectorAffine = MOI.VectorAffineFunction{Float64}

const Interval = MOI.Interval{Float64}
const LessThan = MOI.LessThan{Float64}
const GreaterThan = MOI.GreaterThan{Float64}
const EqualTo = MOI.EqualTo{Float64}
const IntervalConvertible = Union{Interval, LessThan, GreaterThan, EqualTo}

const Zeros = MOI.Zeros
const Nonnegatives = MOI.Nonnegatives
const Nonpositives = MOI.Nonpositives
const SupportedVectorSets = Union{Zeros, Nonnegatives, Nonpositives}

import OSQP

# TODO: consider moving to MOI:
constant(f::MOI.SingleVariable) = 0
constant(f::MOI.ScalarAffineFunction) = f.constant
# constant(f::MOI.ScalarQuadraticFunction) = f.constant

dimension(s::MOI.AbstractSet) = MOI.dimension(s)
dimension(::MOI.AbstractScalarSet) = 1

lower(::Zeros, i::Int) = 0.0
lower(::Nonnegatives, i::Int) = 0.0
lower(::Nonpositives, i::Int) = -Inf
upper(::Zeros, i::Int) = 0.0
upper(::Nonnegatives, i::Int) = Inf
upper(::Nonpositives, i::Int) = 0.0

# TODO: just use ∈ on 0.7 (allocates on 0.6):
function _contains(haystack, needle)
    for x in haystack
        x == needle && return true
    end
    false
end

mutable struct OSQPOptimizer <: MOI.AbstractOptimizer
    inner::OSQP.Model
    hasresults::Bool
    results::OSQP.Results
    isempty::Bool
    settings::Dict{Symbol, Any} # need to store these, because they should be preserved if empty! is called
    sense::MOI.OptimizationSense
    objconstant::Float64
    constrconstant::Vector{Float64}
    modcache::ProblemModificationCache{Float64}
    warmstartcache::WarmStartCache{Float64}
    rowranges::Dict{Int, UnitRange{Int}}

    function OSQPOptimizer()
        inner = OSQP.Model()
        hasresults = false
        results = OSQP.Results()
        isempty = true
        settings = Dict{Symbol, Any}()
        sense = MOI.MinSense
        objconstant = 0.
        constrconstant = Float64[]
        modcache = ProblemModificationCache{Float64}()
        warmstartcache = WarmStartCache{Float64}()
        rowranges = Dict{Int, UnitRange{Int}}()
        new(inner, hasresults, results, isempty, settings, sense, objconstant, constrconstant, modcache, warmstartcache, rowranges)
    end
end

hasresults(optimizer::OSQPOptimizer) = optimizer.hasresults

function MOI.empty!(optimizer::OSQPOptimizer)
    optimizer.inner = OSQP.Model()
    optimizer.hasresults = false
    optimizer.results = OSQP.Results()
    optimizer.isempty = true
    optimizer.sense = MOI.MinSense # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer.constrconstant = Float64[]
    optimizer.modcache = ProblemModificationCache{Float64}()
    optimizer.warmstartcache = WarmStartCache{Float64}()
    empty!(optimizer.rowranges)
    optimizer
end

MOI.isempty(optimizer::OSQPOptimizer) = optimizer.isempty

struct UnsupportedObjectiveError <: Exception end

function MOI.copy!(dest::OSQPOptimizer, src::MOI.ModelLike; copynames=false)
    copynames && error("Copying names is not supported.")
    MOI.empty!(dest)
    idxmap = MOIU.IndexMap(dest, src)
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
    A, l, u, dest.constrconstant = processconstraints(src, idxmap, dest.rowranges)
    OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, dest.settings...)
    dest.modcache = ProblemModificationCache(P, q, A, l, u)
    dest.warmstartcache = WarmStartCache{Float64}(size(A, 2), size(A, 1))
    processprimalstart!(dest.warmstartcache.x, src, idxmap)
    processdualstart!(dest.warmstartcache.y, src, idxmap, dest.rowranges)
    dest.isempty = false
    idxmap
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
        MOI.supportsconstraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F, S}())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = CI{F, S}(i)
        end
    end
    idxmap
end

function assign_constraint_row_ranges!(rowranges::Dict{Int, UnitRange{Int}}, idxmap::MOIU.IndexMap, src::MOI.ModelLike)
    startrow = 1
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci_src in cis_src
            set = MOI.get(src, MOI.ConstraintSet(), ci_src)
            ci_dest = idxmap[ci_src]
            endrow = startrow + dimension(set) - 1
            rowranges[ci_dest.value] = startrow : endrow
            startrow = endrow + 1
        end
    end
end

function constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractScalarSet})
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    first(rowrange)
end
constraint_rows(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractVectorSet}) = rowranges[ci.value]
constraint_rows(optimizer::OSQPOptimizer, ci::CI) = constraint_rows(optimizer.rowranges, ci)

"""
Return objective sense, as well as matrix `P`, vector `q`, and scalar `c` such that objective function is `1/2 x' P x + q' x + c`.
"""
function processobjective(src::MOI.ModelLike, idxmap)
    sense = MOI.get(src, MOI.ObjectiveSense())
    n = MOI.get(src, MOI.NumberOfVariables())
    q = zeros(n)
    if sense != MOI.FeasibilitySense
        if MOI.canget(src, MOI.ObjectiveFunction{MOI.SingleVariable}())
            fsingle = MOI.get(src, MOI.ObjectiveFunction{MOI.SingleVariable}())
            P = spzeros(n, n)
            q[idxmap[fsingle.variable].value] = 1
            c = 0.
        elseif MOI.canget(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
            faffine = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
            P = spzeros(n, n)
            processlinearterms!(q, faffine.terms, idxmap)
            c = faffine.constant
        elseif MOI.canget(src, MOI.ObjectiveFunction{Quadratic}())
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{Quadratic}())
            I = [Int(idxmap[term.variable_index_1].value) for term in fquadratic.quadratic_terms]
            J = [Int(idxmap[term.variable_index_2].value) for term in fquadratic.quadratic_terms]
            V = [term.coefficient for term in fquadratic.quadratic_terms]
            symmetrize!(I, J, V)
            P = sparse(I, J, V, n, n)
            processlinearterms!(q, fquadratic.affine_terms, idxmap)
            c = fquadratic.constant
        else
            throw(UnsupportedObjectiveError())
        end
        sense == MOI.MaxSense && (Compat.rmul!(P, -1); Compat.rmul!(q, -1); c = -c)
    else
        P = spzeros(n, n)
        q = zeros(n)
        c = 0.
    end
    sense, P, q, c
end

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmapfun::Function = identity)
    # This is currently needed to avoid depwarns. TODO: make this nice again:
    if q isa VectorModificationCache
        q[:] = 0
    else
        q .= 0
    end
    for term in terms
        var = term.variable_index
        coeff = term.coefficient
        q[idxmapfun(var).value] += coeff
    end
end

function processlinearterms!(q, terms::Vector{<:MOI.ScalarAffineTerm}, idxmap::MOIU.IndexMap)
    processlinearterms!(q, terms, var -> idxmap[var])
end

function symmetrize!(I::Vector{Int}, J::Vector{Int}, V::Vector)
    n = length(V)
    (length(I) == length(J) == n) || error()
    for i = 1 : n
        if I[i] != J[i]
            push!(I, J[i])
            push!(J, I[i])
            push!(V, V[i])
        end
    end
end

function processconstraints(src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}})
    if VERSION < v"0.7-"
        m = mapreduce(length, +, 0, values(rowranges))
    else
        m = mapreduce(length, +, values(rowranges), init=0)
    end
    l = Vector{Float64}(undef, m)
    u = Vector{Float64}(undef, m)
    constant = Vector{Float64}(undef, m)
    bounds = (l, u)
    I = Int[]
    J = Int[]
    V = Float64[]
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!((I, J, V), bounds, constant, src, idxmap, rowranges, F, S)
    end
    l .-= constant
    u .-= constant
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)
    (A, l, u, constant)
end

function processconstraints!(triplets::SparseTriplets, bounds::Tuple{<:Vector, <:Vector}, constant::Vector{Float64},
        src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_rows(rowranges, idxmap[ci])
        processconstant!(constant, rows, f)
        processlinearpart!(triplets, f, rows, idxmap)
        processconstraintset!(bounds, rows, s)
    end
    nothing
end

function processconstant!(c::Vector{Float64}, row::Int, f::AffineConvertible)
    c[row] = constant(f)
    nothing
end

function processconstant!(c::Vector{Float64}, rows::UnitRange{Int}, f::VectorAffine)
    for (i, row) in enumerate(rows)
        c[row] = f.constants[i]
    end
end

function processlinearpart!(triplets::SparseTriplets, f::MOI.SingleVariable, row::Int, idxmap)
    (I, J, V) = triplets
    col = idxmap[f.variable].value
    push!(I, row)
    push!(J, col)
    push!(V, 1)
    nothing
end

function processlinearpart!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap)
    (I, J, V) = triplets
    for term in f.terms
        var = term.variable_index
        coeff = term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, coeff)
    end
end

function processlinearpart!(triplets::SparseTriplets, f::MOI.VectorAffineFunction, rows::UnitRange{Int}, idxmap)
    (I, J, V) = triplets
    for term in f.terms
        row = rows[term.output_index]
        var = term.scalar_term.variable_index
        coeff = term.scalar_term.coefficient
        col = idxmap[var].value
        push!(I, row)
        push!(J, col)
        push!(V, coeff)
    end
end

function processconstraintset!(bounds::Tuple{<:Vector, <:Vector}, row::Int, s::IntervalConvertible)
    l, u = bounds
    interval = MOI.Interval(s)
    l[row] = interval.lower
    u[row] = interval.upper
    nothing
end

function processconstraintset!(bounds::Tuple{<:Vector, <:Vector}, rows::UnitRange{Int}, s::S) where {S<:SupportedVectorSets}
    l, u = bounds
    for (i, row) in enumerate(rows)
        l[row] = lower(s, i)
        u[row] = upper(s, i)
    end
end

function processprimalstart!(x, src::MOI.ModelLike, idxmap)
    if MOI.canget(src, MOI.VariablePrimalStart(), VI)
        vis_src = MOI.get(src, MOI.ListOfVariableIndices())
        for vi in vis_src
            x[idxmap[vi]] = get(src, MOI.VariablePrimalStart(), vi)
        end
    end
end

function processdualstart!(y, src::MOI.ModelLike, idxmap, rowranges::Dict{Int, UnitRange{Int}})
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if MOI.canget(src, MOI.ConstraintDualStart(), CI{F, S})
            cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            for ci in cis_src
                rows = constraint_rows(rowranges, idxmap[ci])
                dual = get(src, MOI.ConstraintDualStart(), ci)
                for (i, row) in enumerate(rows)
                    y[row] = -dual[i] # opposite dual convention
                end
            end
        end
    end
end


## Standard optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense) = optimizer.sense
MOI.set!(optimizer::OSQPOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense) = optimizer.sense = sense

MOI.canget(optimizer::OSQPOptimizer, ::MOI.NumberOfVariables) = !MOI.isempty(optimizer) # https://github.com/oxfordcontrol/OSQP.jl/issues/10
function MOI.get(optimizer::OSQPOptimizer, a::MOI.NumberOfVariables)
    MOI.canget(optimizer, a) || error()
    OSQP.dimensions(optimizer.inner)[1]
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = MOI.canget(optimizer, MOI.NumberOfVariables())
function MOI.get(optimizer::OSQPOptimizer, a::MOI.ListOfVariableIndices)
    MOI.canget(optimizer, a) || error()
    [VI(i) for i = 1 : MOI.get(optimizer, MOI.NumberOfVariables())] # TODO: support for UnitRange would be nice
end


## Solver-specific optimizer attributes:
module OSQPSettings

using Compat
using MathOptInterface
using OSQP

export OSQPAttribute, isupdatable

abstract type OSQPAttribute <: MathOptInterface.AbstractOptimizerAttribute end

# TODO: just use ∈ on 0.7 (allocates on 0.6):
function _contains(haystack, needle)
    for x in haystack
        x == needle && return true
    end
    false
end

for setting in fieldnames(OSQP.Settings)
    Attribute = Symbol(mapreduce(uppercasefirst, *, split(String(setting), '_'))) # to camelcase
    @eval begin
        export $Attribute
        struct $Attribute <: OSQPAttribute end
        Base.Symbol(::$Attribute) = $(QuoteNode(setting))
        isupdatable(::$Attribute) = $(_contains(OSQP.UPDATABLE_SETTINGS, setting))
    end
end
end # module

using .OSQPSettings

function MOI.set!(optimizer::OSQPOptimizer, a::OSQPAttribute, value)
    (isupdatable(a) || MOI.isempty(optimizer)) || throw(MOI.CannotSetAttribute(a))
    setting = Symbol(a)
    optimizer.settings[setting] = value
    if !MOI.isempty(optimizer)
        OSQP.update_settings!(optimizer.inner; setting => value)
    end
end


## Optimizer methods:
function MOI.optimize!(optimizer::OSQPOptimizer)
    processupdates!(optimizer.inner, optimizer.modcache)
    processupdates!(optimizer.inner, optimizer.warmstartcache)
    OSQP.solve!(optimizer.inner, optimizer.results)
    optimizer.hasresults = true
    # Copy previous solution into warm start cache without setting the dirty bit:
    copyto!(optimizer.warmstartcache.x.data, optimizer.results.x)
    copyto!(optimizer.warmstartcache.y.data, optimizer.results.y)
    nothing
end

# OSQP.Model already sets up a finalizer that calls OSQP.clean!. Manually calling it would result in a double free.
# MOI.free!(optimizer::OSQPOptimizer) = OSQP.clean!(optimizer.inner)


## Optimizer attributes:
MOI.canget(::OSQPOptimizer, ::MOI.RawSolver) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.RawSolver) = optimizer.inner

MOI.canget(::OSQPOptimizer, ::MOI.ResultCount) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ResultCount) = 1

MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{Quadratic}) = true
MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveSense) = true

function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{MOI.SingleVariable}, obj::MOI.SingleVariable)
    MOI.isempty(optimizer) && throw(MOI.CannotSetAttribute(a))
    optimizer.modcache.P[:] = 0
    optimizer.modcache.q[:] = 0
    optimizer.modcache.q[obj.variable.value] = 1
    optimizer.objconstant = 0
    nothing
end

function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}, obj::MOI.ScalarAffineFunction{Float64})
    MOI.isempty(optimizer) && throw(MOI.CannotSetAttribute(a))
    optimizer.modcache.P[:] = 0
    processlinearterms!(optimizer.modcache.q, obj.terms)
    optimizer.objconstant = obj.constant
    nothing
end

function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{Quadratic}, obj::Quadratic)
    MOI.isempty(optimizer) && throw(MOI.CannotSetAttribute(a))
    cache = optimizer.modcache
    cache.P[:] = 0
    for term in obj.quadratic_terms
        row = term.variable_index_1.value
        col = term.variable_index_2.value
        coeff = term.coefficient
        row > col && ((row, col) = (col, row)) # upper triangle only
        cache.P[row, col] += coeff
    end
    processlinearterms!(optimizer.modcache.q, obj.affine_terms)
    optimizer.objconstant = obj.constant
    nothing
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveValue) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, a::MOI.ObjectiveValue)
    MOI.canget(optimizer, a) || error()
    rawobj = optimizer.results.info.obj_val + optimizer.objconstant
    ifelse(optimizer.sense == MOI.MaxSense, -rawobj, rawobj)
end

# Since these aren't explicitly returned by OSQP, I feel like it would be better to have a fallback method compute these:
# MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveBound) = false
# MOI.canget(optimizer::OSQPOptimizer, ::MOI.RelativeGap) = false

MOI.canget(optimizer::OSQPOptimizer, ::MOI.SolveTime) = hasresults(optimizer)
MOI.get(optimizer::OSQPOptimizer, a::MOI.SolveTime) = (MOI.canget(optimizer, a) || error(); optimizer.results.info.run_time)

MOI.canget(optimizer::OSQPOptimizer, ::MOI.TerminationStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, a::MOI.TerminationStatus)
    # Note that the :Dual_infeasible and :Primal_infeasible are mapped to MOI.Success
    # because OSQP can return a proof of infeasibility. For the same reason,
    # :Primal_infeasible_inaccurate is mapped to MOI.AlmostSuccess
    MOI.canget(optimizer, a) || error()
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
    elseif osqpstatus == :Non_convex
        MOI.InvalidModel
    end
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.PrimalStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, a::MOI.PrimalStatus)
    MOI.canget(optimizer, a) || error()
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
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.DualStatus) = hasresults(optimizer)
function MOI.get(optimizer::OSQPOptimizer, a::MOI.DualStatus)
    MOI.canget(optimizer, a) || error()
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
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate, :Non_convex (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end


## Variables:
function MOI.isvalid(optimizer::OSQPOptimizer, vi::VI)
    MOI.canget(optimizer, MOI.NumberOfVariables()) && vi.value ∈ 1 : MOI.get(optimizer, MOI.NumberOfVariables())
end


## Variable attributes:
function MOI.canget(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, ::Type{VI})
    hasresults(optimizer) || return false
    _contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status) || optimizer.results.dual_inf_cert != nothing
end

function MOI.get(optimizer::OSQPOptimizer, a::MOI.VariablePrimal, vi::VI)
    MOI.canget(optimizer, a, typeof(vi)) || error()
    x = ifelse(_contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status), optimizer.results.x, optimizer.results.dual_inf_cert)
    x[vi.value]
end

function MOI.set!(optimizer::OSQPOptimizer, a::MOI.VariablePrimalStart, vi::VI, value)
    MOI.isempty(optimizer) && throw(MOI.CannotSetAttribute(a))
    optimizer.warmstartcache.x[vi.value] = value
end


## Constraints:
function MOI.isvalid(optimizer::OSQPOptimizer, ci::CI)
    MOI.isempty(optimizer) && return false
    ci.value ∈ keys(optimizer.rowranges)
end

function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ConstraintDualStart, ci::CI, value)
    MOI.isempty(optimizer) && throw(MOI.CannotSetAttribute(a))
    rows = constraint_rows(optimizer, ci)
    for (i, row) in enumerate(rows)
        optimizer.warmstartcache.y[row] = -value[i] # opposite dual convention
    end
    nothing
end

# function modification:
function MOI.set!(optimizer::OSQPOptimizer, attr::MOI.ConstraintFunction, ci::CI{Affine, <:IntervalConvertible}, f::Affine)
    MOI.isvalid(optimizer, ci) || error("Invalid constraint index")
    row = constraint_rows(optimizer, ci)
    optimizer.modcache.A[row, :] = 0
    for term in f.terms
        col = term.variable_index.value
        coeff = term.coefficient
        optimizer.modcache.A[row, col] += coeff
    end
    Δconstant = optimizer.constrconstant[row] - f.constant
    optimizer.constrconstant[row] = f.constant
    optimizer.modcache.l[row] += Δconstant
    optimizer.modcache.u[row] += Δconstant
    nothing
end

function MOI.set!(optimizer::OSQPOptimizer, attr::MOI.ConstraintFunction, ci::CI{VectorAffine, <:SupportedVectorSets}, f::VectorAffine)
    MOI.isvalid(optimizer, ci) || error("Invalid constraint index")
    rows = constraint_rows(optimizer, ci)
    for row in rows
        optimizer.modcache.A[row, :] = 0
    end
    for term in f.terms
        row = rows[term.output_index]
        col = term.scalar_term.variable_index.value
        coeff = term.scalar_term.coefficient
        optimizer.modcache.A[row, col] += coeff
    end
    for (i, row) in enumerate(rows)
        Δconstant = optimizer.constrconstant[row] - f.constants[i]
        optimizer.constrconstant[row] = f.constants[i]
        optimizer.modcache.l[row] += Δconstant
        optimizer.modcache.u[row] += Δconstant
    end
end

# set modification:
function MOI.set!(optimizer::OSQPOptimizer, attr::MOI.ConstraintSet, ci::CI{<:AffineConvertible, S}, s::S) where {S <: IntervalConvertible}
    MOI.isvalid(optimizer, ci) || error("Invalid constraint index")
    interval = MOI.Interval(s)
    row = constraint_rows(optimizer, ci)
    constant = optimizer.constrconstant[row]
    optimizer.modcache.l[row] = interval.lower - constant
    optimizer.modcache.u[row] = interval.upper - constant
    nothing
end

function MOI.set!(optimizer::OSQPOptimizer,  attr::MOI.ConstraintSet, ci::CI{<:VectorAffine, S}, s::S) where {S <: SupportedVectorSets}
    MOI.isvalid(optimizer, ci) || error("Invalid constraint index")
    rows = constraint_rows(optimizer, ci)
    for (i, row) in enumerate(rows)
        constant = optimizer.constrconstant[row]
        optimizer.modcache.l[row] = lower(s, i) - constant
        optimizer.modcache.u[row] = upper(s, i) - constant
    end
    nothing
end

# partial function modification:
function MOI.modify!(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, change::MOI.ScalarCoefficientChange)
    MOI.isvalid(optimizer, ci) || error("Invalid constraint index")
    row = constraint_rows(optimizer, ci)
    optimizer.modcache.A[row, change.variable.value] = change.new_coefficient
    nothing
end

# TODO: MultirowChange?

MOI.supportsconstraint(optimizer::OSQPOptimizer, ::Type{<:AffineConvertible}, ::Type{<:IntervalConvertible}) = true
MOI.supportsconstraint(optimizer::OSQPOptimizer, ::Type{VectorAffine}, ::Type{<:SupportedVectorSets}) = true

## Constraint attributes:
function MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDual, ::Type{<:CI})
    hasresults(optimizer) || return false
    _contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status) || optimizer.results.prim_inf_cert != nothing
end

function MOI.get(optimizer::OSQPOptimizer, a::MOI.ConstraintDual, ci::CI)
    MOI.canget(optimizer, a, typeof(ci)) || error()
    y = ifelse(_contains(OSQP.SOLUTION_PRESENT, optimizer.results.info.status), optimizer.results.y, optimizer.results.prim_inf_cert)
    rows = constraint_rows(optimizer, ci)
    -y[rows]
end


# Objective modification
function MOI.modify!(optimizer::OSQPOptimizer, attr::MOI.ObjectiveFunction, change::MOI.ScalarConstantChange)
    MOI.isempty(optimizer) && error()  # TODO: throw a MOI.CannotModifyObjective() exception once that exists
    optimizer.objconstant = change.new_constant
end

function MOI.modify!(optimizer::OSQPOptimizer, attr::MOI.ObjectiveFunction, change::MOI.ScalarCoefficientChange)
    MOI.isempty(optimizer) && error()  # TODO: throw a MOI.CannotModifyObjective() exception once that exists
    optimizer.modcache.q[change.variable.value] = change.new_coefficient
end

# There is currently no ScalarQuadraticCoefficientChange.

MOIU.@model(OSQPModel, # modelname
    (), # scalarsets
    (Interval, LessThan, GreaterThan, EqualTo), # typedscalarsets
    (Zeros, Nonnegatives, Nonpositives), # vectorsets
    (), # typedvectorsets
    (SingleVariable,), # scalarfunctions
    (ScalarAffineFunction, ScalarQuadraticFunction), # typedscalarfunctions
    (), # vectorfunctions
    (VectorAffineFunction,) # typedvectorfunctions
)

end # module
