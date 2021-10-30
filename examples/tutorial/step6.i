[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolids]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3  
    //     E     nu
  [end]
[end]

[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux uy
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dofs=uy
    value=0.02
    boundary=top
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
rank4mate=jacobian
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
  r_rel_tol=5.0e-8
  r_abs_tol=4.5e-7
  solver=mumps
[end]

[job]
  type=static
  debug=dep
[end]