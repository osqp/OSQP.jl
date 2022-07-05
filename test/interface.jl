using SparseArrays: SparseMatrixCSC, sparse
using LinearAlgebra: I
using OSQP: Ccsc, ManagedCcsc
@testset "sparse matrix interface roundtrip" begin
    jl = sparse(Matrix{Bool}(LinearAlgebra.I, 5, 5))
    mc = ManagedCcsc(jl)
    GC.@preserve mc begin
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
