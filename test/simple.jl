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
    println(full(A))
    println(l)
    println(u)

    # Create OSQP Data
    # data = OSQP.OSQPData(n, m, P, q, A, l, u)
    m = OSQP.Model()
    OSQP.setup!(m, P, q, A, l, u)
    results = OSQP.solve(m)


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
