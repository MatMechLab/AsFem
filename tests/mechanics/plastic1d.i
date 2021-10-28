*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=1
  xmax=1
  nx=100
  meshtype=edge2
[end]



[dofs]
name=ux
[end]

[projection]
scalarmate=vonMises effective_plastic_strain
rank2mate=stress strain
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux
    mate=myplastic
    domain=alldomain
  [end]
[end]

[mates]
  [myplastic]
    type=plastic1d
    params=120.0  0.5            1.5
    //     E      yield stress   hardening modulus
  [end]
[end]



[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux
    boundary=left
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dofs=ux
    value=1.0*t
    boundary=right
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-4
  endtime=1.0e-3
  adaptive=false
  optiters=3
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
[end]
