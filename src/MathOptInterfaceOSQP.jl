module MathOptInterfaceOSQP

export OSQPOptimizer, OSQPSettings, OSQPModel

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

mutable struct VectorModificationCache{T}
    data::Vector{T}
    dirty::Bool
    VectorModificationCache(data::Vector{T}) where {T} = new{T}(copy(data), false)
end
Base.setindex!(cache::VectorModificationCache, x, i::Integer) = (cache.dirty = true; cache.data[i] = x)
Base.setindex!(cache::VectorModificationCache, x, ::Colon) = (cache.dirty = true; cache.data[:] = x)
Base.getindex(cache::VectorModificationCache, i::Integer) = cache.data[i]

function processupdates!(model::OSQP.Model, cache::VectorModificationCache, updatefun::Function)
    if cache.dirty
        updatefun(model, cache.data)
        cache.dirty = false
    end
end

struct MatrixModificationCache{T}
    cartesian_indices::Vector{CartesianIndex{2}}
    cartesian_indices_set::Set{CartesianIndex{2}} # to speed up checking whether indices are out of bounds in setindex!
    cartesian_indices_per_row::Dict{Int, Vector{CartesianIndex{2}}}
    modifications::Dict{CartesianIndex{2}, T}
    vals::Vector{T}
    inds::Vector{Int}

    function MatrixModificationCache(S::SparseMatrixCSC{T}) where T
        cartesian_indices = Vector{CartesianIndex{2}}(uninitialized, nnz(S))
        cartesian_indices_per_row = Dict{Int, Vector{CartesianIndex{2}}}()
        sizehint!(cartesian_indices_per_row, nnz(S))
        @inbounds for col = 1 : S.n, k = S.colptr[col] : (S.colptr[col+1]-1) # from sparse findn
            row = S.rowval[k]
            I = CartesianIndex(row, col)
            cartesian_indices[k] = I
            push!(get!(() -> CartesianIndex{2}[], cartesian_indices_per_row, row), I)
        end
        modifications = Dict{CartesianIndex{2}, Int}()
        sizehint!(modifications, nnz(S))
        new{T}(cartesian_indices, Set(cartesian_indices), cartesian_indices_per_row, modifications, T[], Int[])
    end
end

function Base.setindex!(cache::MatrixModificationCache, x, row::Integer, col::Integer)
    I = CartesianIndex(row, col)
    @boundscheck I ∈ cache.cartesian_indices_set || throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    cache.modifications[I] = x
end

function Base.setindex!(cache::MatrixModificationCache, x::Real, ::Colon)
    # used to zero out the entire matrix
    @boundscheck x == 0 || throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    for I in cache.cartesian_indices
        cache.modifications[I] = 0
    end
end

function Base.setindex!(cache::MatrixModificationCache, x::Real, row::Integer, ::Colon)
    # used to zero out a row
    @boundscheck x == 0 || throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    for I in cache.cartesian_indices_per_row[row]
        cache.modifications[I] = 0
    end
end

function Base.getindex(cache::MatrixModificationCache, row::Integer, col::Integer)
    cache.modifications[CartesianIndex(row, col)]
end

function isassigned(cache::MatrixModificationCache, row::Integer, col::Integer)
    haskey(cache.modifications, CartesianIndex(row, col))
end

function processupdates!(model::OSQP.Model, cache::MatrixModificationCache, updatefun::Function)
    dirty = !isempty(cache.modifications)
    if dirty
        nmods = length(cache.modifications)
        resize!(cache.vals, nmods)
        resize!(cache.inds, nmods)
        count = 1
        for i in eachindex(cache.cartesian_indices)
            I = cache.cartesian_indices[i]
            if haskey(cache.modifications, I)
                cache.vals[count] = cache.modifications[I]
                cache.inds[count] = i
                count += 1
            end
        end
        updatefun(model, cache.vals, cache.inds)
        empty!(cache.modifications)
    end
end

struct ProblemModificationCache{T}
    P::MatrixModificationCache{T}
    q::VectorModificationCache{T}
    A::MatrixModificationCache{T}
    l::VectorModificationCache{T}
    u::VectorModificationCache{T}

    ProblemModificationCache{T}() where {T} = new{T}()
    function ProblemModificationCache(P::SparseMatrixCSC, q::Vector{T}, A::SparseMatrixCSC, l::Vector{T}, u::Vector{T}) where T
        MC = MatrixModificationCache
        VC = VectorModificationCache
        new{T}(MC(triu(P)), VC(q), MC(A), VC(l), VC(u))
    end
end

function processupdates!(model::OSQP.Model, cache::ProblemModificationCache)
    processupdates!(model, cache.P, OSQP.update_P!)
    processupdates!(model, cache.q, OSQP.update_q!)
    processupdates!(model, cache.A, OSQP.update_A!)
    processupdates!(model, cache.l, OSQP.update_l!)
    processupdates!(model, cache.u, OSQP.update_u!)
end

struct WarmStartCache{T}
    x::VectorModificationCache{T}
    y::VectorModificationCache{T}

    WarmStartCache{T}() where {T} = new{T}()
    function WarmStartCache{T}(n::Integer, m::Integer) where T
        new{T}(VectorModificationCache(zeros(T, n)), VectorModificationCache(zeros(T, m)))
    end
end

function processupdates!(model::OSQP.Model, cache::WarmStartCache)
    # Ugh. I would greatly prefer just doing this:
    #   processupdates!(model, cache.x, (optimizer, xstart) -> (OSQP.warm_start!(optimizer; x = xstart))) # TODO: non-kwarg function would be preferable
    #   processupdates!(model, cache.y, (optimizer, ystart) -> (OSQP.warm_start!(optimizer; y = ystart))) # TODO: non-kwarg function would be preferable
    # but calling warm_start! with only one of x and y zeros the warm start for the other.
    # For now, just set warm start for both *always*
    OSQP.warm_start!(model; x = cache.x.data, y = cache.y.data)
end

mutable struct OSQPOptimizer <: MOI.AbstractOptimizer
    inner::OSQP.Model
    results::Union{OSQP.Results, Nothing}
    isempty::Bool
    settings::Dict{Symbol, Any} # need to store these, because they should be preserved if empty! is called
    sense::MOI.OptimizationSense
    objconstant::Float64
    modcache::ProblemModificationCache{Float64}
    warmstartcache::WarmStartCache{Float64}

    function OSQPOptimizer()
        new(OSQP.Model(), nothing, true, Dict{Symbol, Any}(), MOI.MinSense, 0., ProblemModificationCache{Float64}(), WarmStartCache{Float64}())
    end
end

hasresults(optimizer::OSQPOptimizer) = optimizer.results != nothing

function MOI.empty!(optimizer::OSQPOptimizer)
    optimizer.inner = OSQP.Model()
    optimizer.results = nothing
    optimizer.isempty = true
    optimizer.sense = MOI.MinSense # model parameter, so needs to be reset
    optimizer.objconstant = 0.
    optimizer.modcache = ProblemModificationCache{Float64}()
    optimizer.warmstartcache = WarmStartCache{Float64}()
    optimizer
end

MOI.isempty(optimizer::OSQPOptimizer) = optimizer.isempty

struct UnsupportedObjectiveError <: Exception end

struct UnsupportedConstraintError
    F::Type
    S::Type
end

function MOI.copy!(dest::OSQPOptimizer, src::MOI.ModelLike)
    try
        MOI.empty!(dest)
        idxmap = MOIU.IndexMap(dest, src)
        dest.sense, P, q, dest.objconstant = processobjective(src, idxmap)
        A, l, u = processconstraints(src, idxmap)
        OSQP.setup!(dest.inner; P = P, q = q, A = A, l = l, u = u, dest.settings...)
        dest.modcache = ProblemModificationCache(P, q, A, l, u)
        dest.warmstartcache = WarmStartCache{Float64}(length(idxmap.varmap), length(idxmap.conmap))
        processprimalstart!(dest.warmstartcache.x, src, idxmap)
        processdualstart!(dest.warmstartcache.y, src, idxmap)
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
            I = [idxmap[var].value for var in fquadratic.quadratic_rowvariables]
            J = [idxmap[var].value for var in fquadratic.quadratic_colvariables]
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
    @assert length(variables) == ncoeffs
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
        processconstraintfun!(triplets, f, row, idxmap)
    end
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

function processdualstart!(y, src::MOI.ModelLike, idxmap)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if MOI.canget(src, MOI.ConstraintDualStart(), CI{F, S})
            cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            for ci in cis_src
                y[idxmap[ci].value] = -get(src, MOI.ConstraintDualStart(), ci) # opposite dual convention
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
    OSQP.dimensions(optimizer.model)[1]
end

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ListOfVariableIndices) = MOI.canget(optimizer, MOI.NumberOfVariables())
function MOI.get(optimizer::OSQPOptimizer, a::MOI.ListOfVariableIndices)
    MOI.canget(optimizer, a) || error()
    [VI(i) for i = 1 : get(optimizer, MOI.NumberOfVariables())] # TODO: support for UnitRange would be nice
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
    optimizer.results = OSQP.solve!(optimizer.inner)
    # TODO: remove this, see comment in processupdates!(model::OSQP.Model, cache::WarmStartCache):
    # Use result as warm start for next solve to mimic what OSQP does normally.
    copy!(optimizer.warmstartcache.x.data, optimizer.results.x)
    copy!(optimizer.warmstartcache.y.data, optimizer.results.y)
end

MOI.free!(optimizer::OSQPOptimizer) = OSQP.clean!(optimizer.inner)


## Optimizer attributes:
MOI.canget(optimizer::OSQPOptimizer, ::MOI.RawSolver) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.RawSolver) = optimizer.inner

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ResultCount) = true
MOI.get(optimizer::OSQPOptimizer, ::MOI.ResultCount) = 1

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
    cache.P[:] = 0
    rows = obj.quadratic_rowvariables
    cols = obj.quadratic_rowvariables
    coeffs = obj.quadratic_coefficients
    n = length(coeffs)
    @assert length(rows) == length(cols) == n

    for i = 1 : n
        row = rows[i].value
        col = cols[i].value
        coeff = coeffs[i]
        row > col && ((row, col) = (col, row)) # upper triangle only
        if isassigned(cache.P, row, col)
            cache.P[row, col] += coeff
        else
            cache.P[row, col] = coeff
        end
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

MOI.canget(optimizer::OSQPOptimizer, ::MOI.ObjectiveBound) = false # TODO
MOI.canget(optimizer::OSQPOptimizer, ::MOI.RelativeGap) = false # TODO

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
MOI.isvalid(optimizer::OSQPOptimizer, vi::VI) = MOI.canget(optimizer, MOI.NumberOfVariables()) && vi.value ∈ 1 : get(optimizer, MOI.NumberOfVariables())
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
    cache = optimizer.warmstartcache
    if hasresults(optimizer) && !cache.x.dirty
        # if warm start hasn't been set yet, update it to last solution to match OSQP's internal behavior
        copy!(cache.x.data, optimizer.results.x)
    end
    cache.x[vi.value] = value
end


## Constraints:
function MOI.isvalid(optimizer::OSQPOptimizer, ci::CI)
    MOI.isempty(optimizer) && return false
    m = OSQP.dimensions(optimizer.inner)[2]
    ci.value ∈ 1 : m
end

MOI.canset(optimizer::OSQPOptimizer, ::MOI.ConstraintDualStart, ::Type{<:CI}) = !MOI.isempty(optimizer)
function MOI.set!(optimizer::OSQPOptimizer, a::MOI.ConstraintDualStart, ci::CI, value)
    MOI.canset(optimizer, a, typeof(ci)) || error()
    cache = optimizer.warmstartcache
    if hasresults(optimizer) && !cache.y.dirty
        # if warm start hasn't been set yet, update it to last solution to match OSQP's internal behavior
        copy!(cache.y.data, optimizer.results.y)
    end
    cache.y[ci.value] = -value # opposite dual convention
end

# function modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, ::Type{Affine}) = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, f::Affine)
    MOI.canmodifyconstraint(optimizer, ci, typeof(f)) || error()
    row = ci.value
    optimizer.modcache.A[row, :] = 0
    ncoeff = length(f.coefficients)
    @assert length(f.variables) == ncoeff
    for i = 1 : ncoeff
        col = f.variables[i].value
        coeff = f.coefficients[i]
        optimizer.modcache.A[row, col] += coeff
    end
end

# set modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{<:AffineConvertible, S}, ::Type{S}) where {S <: IntervalConvertible} = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{<:AffineConvertible, S}, s::S) where {S <: IntervalConvertible}
    MOI.canmodifyconstraint(optimizer, ci, typeof(s)) || error()
    interval = MOI.Interval(s)
    row = ci.value
    optimizer.modcache.l[row] = interval.lower
    optimizer.modcache.u[row] = interval.upper
end

# partial function modification:
MOI.canmodifyconstraint(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, ::Type{MOI.ScalarCoefficientChange{<:Real}}) = MOI.isvalid(optimizer, ci)
function MOI.modifyconstraint!(optimizer::OSQPOptimizer, ci::CI{Affine, <:IntervalConvertible}, change::MOI.ScalarCoefficientChange{<:Real})
    MOI.canmodifyconstraint(optimizer, ci, typeof(change)) || error()
    optimizer.modcache.A[ci.value, change.variable.value] = change.new_coefficient
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
    -y[ci.value] # MOI uses opposite dual convention
end


# Objective modification
MOI.canmodifyobjective(optimizer::OSQPOptimizer, ::Type{<:MOI.ScalarConstantChange}) = !MOI.isempty(optimizer)
function MOI.modifyobjective!(optimizer::OSQPOptimizer, change::MOI.ScalarConstantChange)
    MOI.canmodifyobjective(optimizer, typeof(change)) || error()
    optimizer.objconst = change.new_constant
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
