# OSQP.jl

[![Build Status](https://github.com/oxfordcontrol/OSQP.jl/workflows/CI/badge.svg)](https://github.com/oxfordcontrol/OSQP.jl/actions) 
[![codecov.io](http://codecov.io/github/oxfordcontrol/OSQP.jl/coverage.svg?branch=master)](http://codecov.io/github/oxfordcontrol/OSQP.jl?branch=master)

Julia wrapper for [OSQP](https://osqp.org/): the Operator Splitting QP Solver.

The OSQP (Operator Splitting Quadratic Program) solver is a numerical optimization package for solving problems in the form
```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```

where `x in R^n` is the optimization variable. The objective function is defined by a positive semidefinite matrix `P in S^n_+` and vector `q in R^n`. The linear constraints are defined by matrix `A in R^{m x n}` and vectors `l in R^m U {-inf}^m`, `u in R^m U {+inf}^m`.


## Documentation
The interface is documented [here](https://osqp.org/).
