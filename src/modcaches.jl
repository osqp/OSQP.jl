module ModificationCaches

import LinearAlgebra
import OSQP
import SparseArrays

export VectorModificationCache,
    MatrixModificationCache,
    ProblemModificationCache,
    WarmStartCache,
    processupdates!

mutable struct VectorModificationCache{T}
    data::Vector{T}
    dirty::Bool
    function VectorModificationCache(data::Vector{T}) where {T}
        return new{T}(copy(data), false)
    end
end
function Base.setindex!(cache::VectorModificationCache, x, i::Integer)
    cache.dirty = true
    cache.data[i] = x
    return x
end

function Base.setindex!(cache::VectorModificationCache, x, ::Colon)
    cache.dirty = true
    cache.data .= x
    return x
end

Base.getindex(cache::VectorModificationCache, i::Integer) = cache.data[i]

function processupdates!(
    model,
    cache::VectorModificationCache,
    updatefun::Function,
)
    if cache.dirty
        updatefun(model, cache.data)
        cache.dirty = false
    end
    return
end

struct MatrixModificationCache{T}
    cartesian_indices::Vector{CartesianIndex{2}}
    # To speed up checking whether indices are out of bounds in setindex!
    cartesian_indices_set::Set{CartesianIndex{2}}
    cartesian_indices_per_row::Dict{Int,Vector{CartesianIndex{2}}}
    modifications::Dict{CartesianIndex{2},T}
    vals::Vector{T}
    inds::Vector{Int}

    function MatrixModificationCache(
        S::SparseArrays.SparseMatrixCSC{T},
    ) where {T}
        cartesian_indices =
            Vector{CartesianIndex{2}}(undef, SparseArrays.nnz(S))
        cartesian_indices_per_row = Dict{Int,Vector{CartesianIndex{2}}}()
        sizehint!(cartesian_indices_per_row, SparseArrays.nnz(S))
        @inbounds for col in 1:S.n
            for k in S.colptr[col]:(S.colptr[col+1]-1) # from sparse findn
                row = S.rowval[k]
                I = CartesianIndex(row, col)
                cartesian_indices[k] = I
                if !haskey(cartesian_indices_per_row, row)
                    cartesian_indices_per_row[row] = CartesianIndex{2}[]
                end
                push!(cartesian_indices_per_row[row], I)
            end
        end
        modifications = Dict{CartesianIndex{2},Int}()
        sizehint!(modifications, SparseArrays.nnz(S))
        return new{T}(
            cartesian_indices,
            Set(cartesian_indices),
            cartesian_indices_per_row,
            modifications,
            T[],
            Int[],
        )
    end
end

function Base.setindex!(
    cache::MatrixModificationCache,
    x,
    row::Integer,
    col::Integer,
)
    I = CartesianIndex(row, col)
    @boundscheck if !(I ∈ cache.cartesian_indices_set)
        throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    end
    cache.modifications[I] = x
    return x
end

function Base.setindex!(cache::MatrixModificationCache, x::Real, ::Colon)
    # used to zero out the entire matrix
    @boundscheck if x != 0
        throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    end
    for I in cache.cartesian_indices
        cache.modifications[I] = 0
    end
    return x
end

function Base.setindex!(
    cache::MatrixModificationCache,
    x::Real,
    row::Integer,
    ::Colon,
)
    # used to zero out a row
    @boundscheck if x != 0
        throw(ArgumentError("Changing the sparsity pattern is not allowed."))
    end
    for I in cache.cartesian_indices_per_row[row]
        cache.modifications[I] = 0
    end
    return x
end

function Base.getindex(
    cache::MatrixModificationCache,
    row::Integer,
    col::Integer,
)
    return cache.modifications[CartesianIndex(row, col)]
end

function processupdates!(
    model,
    cache::MatrixModificationCache,
    updatefun::Function,
)
    if !isempty(cache.modifications)
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
    return
end

# More OSQP-specific from here on:
struct ProblemModificationCache{T}
    P::MatrixModificationCache{T}
    q::VectorModificationCache{T}
    A::MatrixModificationCache{T}
    l::VectorModificationCache{T}
    u::VectorModificationCache{T}

    ProblemModificationCache{T}() where {T} = new{T}()
    function ProblemModificationCache(
        P::SparseArrays.SparseMatrixCSC,
        q::Vector{T},
        A::SparseArrays.SparseMatrixCSC,
        l::Vector{T},
        u::Vector{T},
    ) where {T}
        MC = MatrixModificationCache
        VC = VectorModificationCache
        return new{T}(MC(LinearAlgebra.triu(P)), VC(q), MC(A), VC(l), VC(u))
    end
end

function processupdates!(model::OSQP.Model, cache::ProblemModificationCache)
    if cache.l.dirty && cache.u.dirty
        # Special case because setting just l or u may cause an 'upper bound
        # must be greater than or equal to lower bound' error
        OSQP.update_bounds!(model, cache.l.data, cache.u.data)
        cache.l.dirty = false
        cache.u.dirty = false
    end
    processupdates!(model, cache.P, OSQP.update_P!)
    processupdates!(model, cache.q, OSQP.update_q!)
    processupdates!(model, cache.A, OSQP.update_A!)
    processupdates!(model, cache.l, OSQP.update_l!)
    processupdates!(model, cache.u, OSQP.update_u!)
    return
end

struct WarmStartCache{T}
    x::VectorModificationCache{T}
    y::VectorModificationCache{T}

    WarmStartCache{T}() where {T} = new{T}()
    function WarmStartCache{T}(n::Integer, m::Integer) where {T}
        return new{T}(
            VectorModificationCache(zeros(T, n)),
            VectorModificationCache(zeros(T, m)),
        )
    end
end

function processupdates!(model::OSQP.Model, cache::WarmStartCache)
    if cache.x.dirty && cache.y.dirty
        # Special case because setting warm start for x only zeroes the stored
        # warm start for y and vice versa.
        OSQP.warm_start_x_y!(model, cache.x.data, cache.y.data)
        cache.x.dirty = false
        cache.y.dirty = false
    end
    processupdates!(model, cache.x, OSQP.warm_start_x!)
    processupdates!(model, cache.y, OSQP.warm_start_y!)
    return
end

end
