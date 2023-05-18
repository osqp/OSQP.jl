# OSQP.jl

[![Build Status](https://github.com/osqp/OSQP.jl/workflows/CI/badge.svg)](https://github.com/osqp/OSQP.jl/actions) 
[![codecov.io](http://codecov.io/github/osqp/OSQP.jl/coverage.svg?branch=master)](http://codecov.io/github/osqp/OSQP.jl?branch=master)

[OSQP.jl](https://github.com/osqp/OSQP.jl) is a Julia wrapper for
[OSQP](https://osqp.org/): the Operator Splitting QP Solver.

## License

OSQP.jl is licensed under the [Apache-2.0 license](https://github.com/osqp/OSQP.jl/blob/master/LICENSE.md).

The upstream solver, [osqp/osqp](https://github.com/osqp/osqp) is also licensed
under the [Apache-2.0 license](https://github.com/osqp/osqp/blob/master/LICENSE).

## Installation

Install OSQP.jl using the Julia package manager

```julia
import Pkg
Pkg.add("OSQP")
```

## Problem class

The OSQP (Operator Splitting Quadratic Program) solver is a numerical
optimization package for solving problems in the form
```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```
where `x in R^n` is the optimization variable. The objective function is defined
by a positive semidefinite matrix `P in S^n_+` and vector `q in R^n`. The linear
constraints are defined by matrix `A in R^{m x n}` and vectors
`l in R^m U {-inf}^m`, `u in R^m U {+inf}^m`.

## Documentation

Detailed documentation is available at [https://osqp.org/](https://osqp.org/).
