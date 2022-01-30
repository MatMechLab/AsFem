*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=2
  xmax=5.0
  ymax=5.0
  nx=50
  ny=50
  meshtype=quad4
[end]



[dofs]
name=ux uy
[end]

[projection]
scalarmate=vonMises effective_plastic_strain
rank2mate=stress strain
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux uy
    mate=myplastic
    domain=alldomain
  [end]
[end]

[mates]
  [myplastic]
    type=j2plasticity
    params=210.0  0.3  0.5            1.2
    //     E      nu   yield stress   hardening modulus
  [end]
[end]



[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux
    boundary=left right
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dofs=uy
    boundary=bottom
    value=0.0
  [end]
  [loadUy]
    type=cyclicdirichlet
    dofs=uy
    params=0.0 0.0 0.2 0.2 0.3 0.0 0.4 0.4 0.5 -0.1 0.6 0.5
    boundary=top
  [end]
[end]

[output]
type=vtu
interval=2
[end]

[timestepping]
  type=be
  dt=2.0e-4
  endtime=1.0e0
  adaptive=true
  optiters=3
  dtmax=1.0e-1
  //dtmin=1.0e-12
[end]

[job]
  type=transient
  debug=dep
[end]
