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
    dof=ux
    boundary=left right
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

[output]
type=vtu
interval=5
[end]

[timestepping]
  type=be
  dt=2.0e-4
  endtime=3.0e-4
  adaptive=true
  optiters=3
  dtmax=1.0e-1
  dtmin=1.0e-4
[end]

[job]
  type=transient
  debug=dep
[end]
