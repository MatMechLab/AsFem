[mesh]
  type=asfem
  dim=2
  xmax=10.0
  ymax=2.0
  nx=100
  ny=20
  meshtype=quad9
[end]

[dofs]
name=T
[end]

[elmts]
  [mysolids]
    type=thermal
    dofs=T
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=thermalmate
    params=1.0 1.5 2.0 0.0
    //     rho Cp  K   Q
  [end]
[end]

[bcs]
  [flux]
    type=neumann
    dofs=T
    value=-0.1
    boundary=right
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=1.0e-2
  time=1.0e-1
  optiters=3
  adaptive=true
[end]

[projection]
vectormate=gradT
[end]

[job]
  type=transient
  debug=dep
[end]

