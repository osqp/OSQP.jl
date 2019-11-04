using SparseArrays: SparseMatrixCSC, sparse
using LinearAlgebra: I
using OSQP: Ccsc, ManagedCcsc
@testset "sparse matrix interface roundtrip" begin
    jl = sparse(I(5))
    mc = ManagedCcsc(jl)
    c = Ccsc(mc)
    jl2 = Base.unsafe_convert(SparseMatrixCSC, c)
    @test jl == jl2
end
