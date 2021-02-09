using SparseArrays: SparseMatrixCSC, sparse
using LinearAlgebra: I
using OSQP: Ccsc, ManagedCcsc

if Sys.WORD_SIZE != 32 # FIXME fails for 32 bits
@testset "sparse matrix interface roundtrip" begin
    jl = sparse(Matrix{Bool}(LinearAlgebra.I, 5, 5))
    mc = ManagedCcsc(jl)
    c = Ccsc(mc)
    jl2 = convert(SparseMatrixCSC, c)
    @test jl == jl2
end
end

# Check model error handling
@testset "Model error handling" begin
    model = OSQP.Model()
    @test_throws ErrorException OSQP.solve!(model)
end
