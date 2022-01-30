[mesh]
  type=asfem
  dim=3
  xmax=2.0
  ymax=2.0
  zmax=10.0
  nx=5
  ny=5
  nz=50
  meshtype=hex8
  savemesh=true
[end]

[qpoints]
  type=gauss
  order=4
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
    boundary=front
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

