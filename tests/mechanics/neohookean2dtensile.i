*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=2
  xmax=2
  ymax=2
  nx=80
  ny=80
  meshtype=quad4
[end]



[dofs]
name=ux uy
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
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
    boundary=left
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-3
  endtime=2.0e-3
  adaptive=false
  optiters=3
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
[end]
