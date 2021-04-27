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
    mate=smallstrain
    domain=alldomain
  [end]
[end]

[mates]
  [smallstrain]
    type=incresmallstrain
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
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=ux
    value=1.0*t
    boundary=top
  [end]
  [loadUy]
    type=dirichlet
    dof=uy
    value=-2.0*t
    boundary=top
  [end]
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
