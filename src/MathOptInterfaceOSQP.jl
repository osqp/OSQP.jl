module MathOptInterfaceOSQP

export OSQPOptimizer, OSQPSettings, OSQPModel

include("modcaches.jl")
using .ModificationCaches

using Compat
using MathOptInterface
using MathOptInterface.Utilities

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const SparseTriplets = Tuple{Vector{Int}, Vector{Int}, Vector{<:Any}}

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
# constant(f::Quadratic) = f.constant

dimension(s::MOI.AbstractSet) = MOI.dimension(s)
dimension(::MOI.AbstractScalarSet) = 1

mutable struct OSQPOptimizer <: MOI.AbstractOptimizer
    inner::OSQP.Model
    hasresults::Bool
    results::OSQP.Results
    isempty::Bool
    settings::Dict{Symbol, Any} # need to store these, because they should be preserved if empty! is called
    sense::MOI.OptimizationSense
    objconstant::Float64
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
        modcache = ProblemModificationCache{Float64}()
        warmstartcache = WarmStartCache{Float64}()
        rowranges = Dict{Int, UnitRange{Int}}()
        new(inner, hasresults, results, isempty, settings, sense, objconstant, modcache, warmstartcache, rowranges)
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
    optimizer.modcache = ProblemModificationCache{Float64}()
    optimizer.warmstartcache = WarmStartCache{Float64}()
    empty!(optimizer.rowranges)
    optimizer
end

MOI.isempty(optimizer::OSQPOptimizer) = optimizer.isempty

struct UnsupportedObjectiveError <: Exception end

struct UnsupportedConstraintError
    F::Type
    S::Type
end

function MOI.copy!(dest::OSQPOptimizer, src::MOI.ModelLike; copynames=false)
    copynames && error("Copying names is not supported.")
    try
        MOI.empty!(dest)
        idxmap = MOIU.IndexMap(dest, src)
        assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
        dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
        A, l, u = processconstraints(src, idxmap, dest.rowranges)
        OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, dest.settings...)
        dest.modcache = ProblemModificationCache(P, q, A, l, u)
        dest.warmstartcache = WarmStartCache{Float64}(length(idxmap.varmap), length(idxmap.conmap))
        processprimalstart!(dest.warmstartcache.x, src, idxmap)
        processdualstart!(dest.warmstartcache.y, src, idxmap, dest.rowranges)
        dest.isempty = false
        return MOI.CopyResult(MOI.CopySuccess, "", idxmap)
    catch e
        e isa UnsupportedObjectiveError && return MOI.CopyResult(MOI.CopyOtherError, "Unsupported objective", MOIU.IndexMap())
        e isa UnsupportedConstraintError && return MOI.CopyResult(MOI.CopyUnsupportedConstraint, "Unsupported $(e.F)-in-$(e.S) constraint", MOIU.IndexMap())
        rethrow(e)
    end
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
        MOI.supportsconstraint(dest, F, S) || throw(UnsupportedConstraintError(F, S))
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

function constraint_row_indices(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractScalarSet})
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    first(rowrange)
end
constraint_row_indices(rowranges::Dict{Int, UnitRange{Int}}, ci::CI{<:Any, <:MOI.AbstractVectorSet}) = rowranges[ci.value]
constraint_row_indices(optimizer::OSQPOptimizer, ci::CI) = constraint_row_indices(optimizer.rowranges, ci)

"""
Return objective sense, as well as matrix `P`, vector `q`, and scalar `c` such that objective function is 1/2 `x' P x + q' x + c`.
"""
function processobjective(src::MOI.ModelLike, idxmap)
    sense = MOI.get(src, MOI.ObjectiveSense())
    n = MOI.get(src, MOI.NumberOfVariables())
    q = zeros(n)
    if sense != MOI.FeasibilitySense
        if MOI.canget(src, MOI.ObjectiveFunction{SingleVariable}())
            fsingle = MOI.get(src, MOI.ObjectiveFunction{SingleVariable}())
            P = spzeros(n, n)
            q[idxmap[fsingle.variable].value] = 1
            c = 0.
        elseif MOI.canget(src, MOI.ObjectiveFunction{Affine}())
            faffine = MOI.get(src, MOI.ObjectiveFunction{Affine}())
            P = spzeros(n, n)
            processlinearterms!(q, faffine.variables, faffine.coefficients, idxmap)
            c = faffine.constant
        elseif MOI.canget(src, MOI.ObjectiveFunction{Quadratic}())
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{Quadratic}())
            I = [Int(idxmap[var].value) for var in fquadratic.quadratic_rowvariables]
            J = [Int(idxmap[var].value) for var in fquadratic.quadratic_colvariables]
            V = fquadratic.quadratic_coefficients
            symmetrize!(I, J, V)
            P = sparse(I, J, V, n, n)
            processlinearterms!(q, fquadratic.affine_variables, fquadratic.affine_coefficients, idxmap)
            c = fquadratic.constant
        else
            throw(UnsupportedObjectiveError())
        end
        sense == MOI.MaxSense && (scale!(P, -1); scale!(q, -1); c = -c)
    else
        P = spzeros(n, n)
        q = zeros(n)
        c = 0.
    end
    sense, P, q, c
end

function processlinearterms!(q, variables::Vector{VI}, coefficients::Vector, idxmapfun::Function = identity)
    q[:] = 0
    ncoeffs = length(coefficients)
    length(variables) == ncoeffs || error()
    for i = 1 : ncoeffs
        var = variables[i]
        coeff = coefficients[i]
        q[idxmapfun(var).value] += coeff
    end
end

function processlinearterms!(q, variables::Vector{VI}, coefficients::Vector, idxmap::MOIU.IndexMap)
    processlinearterms!(q, variables, coefficients, var -> idxmap[var])
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
    m = reduce(+, 0, length(range) for range in values(rowranges))
    l = Vector{Float64}(undef, m)
    u = Vector{Float64}(undef, m)
    bounds = (l, u)

    I = Int[]
    J = Int[]
    V = Float64[]
    triplets = (I, J, V)

    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        processconstraints!(triplets, bounds, src, idxmap, rowranges, F, S)
    end
    n = MOI.get(src, MOI.NumberOfVariables())
    A = sparse(I, J, V, m, n)

    (A, l, u)
end

function processconstraints!(triplets::SparseTriplets, bounds::Tuple{<:Vector, <:Vector}, src::MOI.ModelLike,
        idxmap, rowranges::Dict{Int, UnitRange{Int}},
        F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_row_indices(rowranges, idxmap[ci])
        processconstraintbounds!(bounds, rows, f, s)
        processconstraintfun!(triplets, f, rows, idxmap)
    end
end

function processconstraintbounds!(bounds::Tuple{<:Vector, <:Vector}, row::Int, f::MOI.AbstractFunction, s::MOI.AbstractScalarSet)
    l, u = bounds
    interval = MOI.Interval(s)
    l[row] = interval.lower - constant(f)
    u[row] = interval.upper - constant(f)
    nothing
end

function processconstraintfun!(triplets::SparseTriplets, f::MOI.SingleVariable, row::Int, idxmap)
    (I, J, V) = triplets
    col = idxmap[f.variable].value
    push!(I, row)
    push!(J, col)
    push!(V, 1)
    nothing
end

function processconstraintfun!(triplets::SparseTriplets, f::MOI.ScalarAffineFunction, row::Int, idxmap)
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
                rows = constraint_row_indices(rowranges, idxmap[ci])
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

export OSQPAttribute, isupdatable

using MathOptInterface, OSQP

abstract type OSQPAttribute <: MathOptInterface.AbstractOptimizerAttribute end

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

MOI.canset(optimizer::OSQPOptimizer, a::OSQPAttribute) = isupdatable(a) || MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::OSQPAttribute, value)
    MOI.canset(optimizer, a) || error()
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
    copy!(optimizer.warmstartcache.x.data, optimizer.results.x)
    copy!(optimizer.warmstartcache.y.data, optimizer.results.y)
    nothing
end

# OSQP.Model already sets up a finalizer that calls OSQP.clean!. Manually calling it would result in a double free.
# MOI.free!(optimizer::OSQPOptimizer) = OSQP.clean!(optimizer.inner)


## Optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.RawSolver) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.RawSolver) = optimizer.inner

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ResultCount) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ResultCount) = 1

MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{SingleVariable}) = true
MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{Affine}) = true
MOI.supports(::OSQPOptimizer, ::MOI.ObjectiveFunction{Quadratic}) = true

MOI.canset(optimizer::OSQPOptimizer, ::MOI.ObjectiveFunction{SingleVariable}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{SingleVariable}, obj::SingleVariable)
    MOI.canset(optimizer, a) || error()
    optimizer.modcache.P[:] = 0
    optimizer.modcache.q[:] = 0
    optimizer.modcache.q[obj.variable.value] = 1
    optimizer.objconstant = 0
    nothing
end

MOI.canset(optimizer::OSQPOptimizer, ::MOI.ObjectiveFunction{Affine}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{Affine}, obj::Affine)
    MOI.canset(optimizer, a) || error()
    optimizer.modcache.P[:] = 0
    processlinearterms!(optimizer.modcache.q, obj.variables, obj.coefficients)
    optimizer.objconstant = obj.constant
    nothing
end

MOI.canset(optimizer::OSQPOptimizer, ::MOI.ObjectiveFunction{Quadratic}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ObjectiveFunction{Quadratic}, obj::Quadratic)
    MOI.canset(optimizer, a) || error()
    cache = optimizer.modcache
    rows = obj.quadratic_rowvariables
    cols = obj.quadratic_rowvariables
    coeffs = obj.quadratic_coefficients
    n = length(coeffs)
    @assert length(rows) == length(cols) == n

    cache.P[:] = 0
    for i = 1 : n
        row = rows[i].value
        col = cols[i].value
        coeff = coeffs[i]
        row > col && ((row, col) = (col, row)) # upper triangle only
        cache.P[row, col] += coeff
    end
    processlinearterms!(optimizer.modcache.q, obj.affine_variables, obj.affine_coefficients)
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
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
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
    else # :Interrupted, :Max_iter_reached, :Solved_inaccurate (TODO: good idea? use OSQP.SOLUTION_PRESENT?)
        MOI.UnknownResultStatus
    end
end


## Variables:
MOI.isvalid(optimizer::OSQPOptimizer, vi::VI) = MOI.canget(optimizer, MOI.NumberOfVariables()) && vi.value ∈ 1 : MOI.get(optimizer, MOI.NumberOfVariables())
MOI.canaddvariable(optimizer::OSQPOptimizer) = false


## Variable attributes:
function MOI.canget(optimizer::OSQPOptimizer, ::MOI.VariablePrimal, ::Type{VI})
    hasresults(optimizer) || return false
    optimizer.results.info.status ∈ OSQP.SOLUTION_PRESENT || optimizer.results.dual_inf_cert != nothing
end

function MOI.get(optimizer::OSQPOptimizer, a::MOI.VariablePrimal, vi::VI)
    MOI.canget(optimizer, a, typeof(vi)) || error()
    x = ifelse(optimizer.results.info.status ∈ OSQP.SOLUTION_PRESENT, optimizer.results.x, optimizer.results.dual_inf_cert)
    x[vi.value]
end

MOI.canset(optimizer::OSQPOptimizer, ::MOI.VariablePrimalStart, ::Type{VI}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.VariablePrimalStart, vi::VI, value)
    MOI.canset(optimizer, a, typeof(vi)) || error()
    optimizer.warmstartcache.x[vi.value] = value
end


## Constraints:
function MOI.isvalid(optimizer::OSQPOptimizer, ci::CI)
    MOI.isempty(optimizer) && return false
    ci.value ∈ keys(optimizer.rowranges)
end

MOI.canset(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{<:CI}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ConstraintDualStart, ci::CI, value)
    MOI.canset(optimizer, a, typeof(ci)) || error()
    rows = constraint_row_indices(optimizer, ci)
    for (i, row) in enumerate(rows)
        optimizer.warmstartcache.y[row] = -value[i]
    end
    nothing
end

# function modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, ::Type{Affine}) = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, f::Affine) # TODO: generalize to vector constraints
    MOI.canmodifyconstraint(optimizer, ci, typeof(f)) || error()
    ncoeff = length(f.coefficients)
    length(f.variables) == ncoeff || error()
    rowrange = constraint_row_indices(optimizer, ci)
    for row in rowrange
        optimizer.modcache.A[row, :] = 0
        for i = 1 : ncoeff
            col = f.variables[i].value
            coeff = f.coefficients[i]
            optimizer.modcache.A[row, col] += coeff
        end
    end
end

# set modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{<:AffineConvertible, S}, ::Type{S}) where {S <: IntervalConvertible} = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{<:AffineConvertible, S}, s::S) where {S <: IntervalConvertible}
    MOI.canmodifyconstraint(optimizer, ci, typeof(s)) || error()
    interval = MOI.Interval(s)
    row = constraint_row_indices(optimizer, ci)
    optimizer.modcache.l[row] = interval.lower
    optimizer.modcache.u[row] = interval.upper
    nothing
end

# partial function modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, ::Type{<:MOI.ScalarCoefficientChange}) = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, change::MOI.ScalarCoefficientChange)
    MOI.canmodifyconstraint(optimizer, ci, typeof(change)) || error()
    row = constraint_row_indices(optimizer, ci)
    optimizer.modcache.A[row, change.variable.value] = change.new_coefficient
    nothing
end

# TODO: MultirowChange?

MOI.supportsconstraint(optimizer::OSQPOptimizer, ::Type{<:AffineConvertible}, ::Type{<:IntervalConvertible}) = true


## Constraint attributes:
function MOI.canget(optimizer::OSQPOptimizer, ::MOI.ConstraintDual, ::Type{<:CI})
    hasresults(optimizer) || return false
    optimizer.results.info.status ∈ OSQP.SOLUTION_PRESENT || optimizer.results.prim_inf_cert != nothing
end

function MOI.get(optimizer::OSQPOptimizer, a::MOI.ConstraintDual, ci::CI)
    MOI.canget(optimizer, a, typeof(ci)) || error()
    y = ifelse(optimizer.results.info.status ∈ OSQP.SOLUTION_PRESENT, optimizer.results.y, optimizer.results.prim_inf_cert)
    rows = constraint_row_indices(optimizer, ci)
    -y[rows]
end


# Objective modification
MOI.canmodifyobjective(optimizer::OSQPOptimizer, ::Type{<:MOI.ScalarConstantChange}) = !MOI.isempty(optimizer)
function MOI.modifyobjective!(optimizer::OSQPOptimizer, change::MOI.ScalarConstantChange)
    MOI.canmodifyobjective(optimizer, typeof(change)) || error()
    optimizer.objconstant = change.new_constant
end

MOI.canmodifyobjective(optimizer::OSQPOptimizer, ::Type{<:MOI.ScalarCoefficientChange}) = !MOI.isempty(optimizer)
function MOI.modifyobjective!(optimizer::OSQPOptimizer, change::MOI.ScalarCoefficientChange)
    MOI.canmodifyobjective(optimizer, typeof(change)) || error()
    optimizer.modcache.q[change.variable.value] = change.new_coefficient
end

# There is currently no ScalarQuadraticCoefficientChange.

MOIU.@model(OSQPModel, # modelname
    (), # scalarsets
    (Interval, LessThan, GreaterThan, EqualTo), # typedscalarsets
    (), # vectorsets
    (), # typedvectorsets
    (SingleVariable,), # scalarfunctions
    (ScalarAffineFunction, ScalarQuadraticFunction), # typedscalarfunctions
    (), # vectorfunctions
    () # typedvectorfunctions
)

end # module
