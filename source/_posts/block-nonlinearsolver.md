---
title: nonlinear block
date: 2021-01-16 21:28:08
categories:
- Document
tags:
- blocks
- input file
- nonlinear solver
---

# [nonlinearsolver] block
The block `[nonlinearsolver]` is used for the newton-raphson nonlinear iteration. This block's layout looks as follows:
```
[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-8
  r_abs_tol=1.0e-7
  stol=1.0e-16
  solver=superlu
[end]
```
## options
The `type =` option specifies the name of the nonlinear solver method one wants to use.

`maxiters=` specifies maximum nonlinear iterations.

`r_rel_tol` and `r_abs_tol` represent the relative and absolute errors of the residual($||R||$).

`stol=` determines tolerance of the displacement error, namely $||\Delta U||$.

`solver=` indicates the name of linear solver, it could be: `ksp`, `mumps`, and `superlu`.


### supported nonlinear solver method
The full list of the available nonlinear solver method is:
```
type=nr // for newton-raphson with line search, the default one
type=newtonls // same as the first one
type=newtontr // NR with trust region method
type=bfgs     // BFGS method
type=broyden  // Broyden method
```
