*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=2
  xmax=5
  ymax=5
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain PK1 PK2
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux uy
    mate=neohookean
    domain=alldomain
  [end]
[end]

[mates]
  [neohookean]
    type=neohookean
    params=100.0 0.3
  [end]
[end]

[nonlinearsolver]
  type=nr
  r_abs_tol=1.0e-10
  r_rel_tol=5.0e-15
  //solver=mumps
[end]

[bcs]
  [FixUx]
    type=dirichlet
    dof=ux
    boundary=bottom
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dof=uy
    boundary=bottom top
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=ux
    value=2.0e-2
    boundary=top
  [end]
[end]



[job]
  type=static
  debug=dep
[end]
