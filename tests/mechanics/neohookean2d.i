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
rank2mate=stress strain PK1 PK2 F
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
    value=1.0*t
    boundary=top
  [end]
[end]

[nonlinearsolver]
  type=nr
  r_abs_tol=4.0e-7
  r_rel_tol=5.0e-8
  //solver=mumps
[end]

[timestepping]
  type=be
  dt=1.0e-3
  endtime=5.0e-3
  adaptive=true
  optiters=3
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
[end]
