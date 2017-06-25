using OSQP

function general_tests()

    # Dimensions
    n = 2
    m = 3
    density = 0.6

    # Cost 
    P = speye(n)
    q = randn(n)
 
    # Constraints
    A = sprandn(m, n, density) 
    l = randn(m)
    u = l + 2 * rand(m)


    # DEBUG
    # println(full(A))
    # println(l)
    # println(u)

    # Create OSQP Data
    # data = OSQP.OSQPData(n, m, P, q, A, l, u)
    model = OSQP.Model()
    # OSQP.setup!(m, P, q, A, l, u; verbose=false)
    OSQP.setup!(model, P, q, A, l, u, verbose=true, alpha=1.8)
    results = OSQP.solve!(model)


    # Update linear cost
    q_new = randn(n) * 2
    OSQP.update!(model, q=q_new)
    results_new = OSQP.solve!(model)


    # Update Px elements
    nnzP = length(P.nzval)
    Px = 10 * randn(nnzP)
    OSQP.update!(model, Px=Px)
    results = OSQP.solve!(model) 



end

out = general_tests()

# function test_basic_qp()
#
#     # Define problem
#
#     # Setup
#
#     # Solve
#
#     # Cleanup
#
#
# end
