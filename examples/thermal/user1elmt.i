[mesh]
  type=asfem
  dim=2
  xmax=5
  ymax=1
  nx=100
  ny=20
  meshtype=quad9
[end]

[dofs]
name=T
[end]

[elmts]
  [mythermal]
    type=user1
    dofs=T
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=user2
    params=1.0 2.0 0.05 1.0
    //     rho Cp  K    Q
  [end]
[end]

[bcs]
  [flux]
    type=neumann
    dofs=T
    value=-0.01
    boundary=right
  [end]
[end]


[timestepping]
  type=be
  dt=1.0e-5
  dtmax=0.1
  time=1.0e-1
  adaptive=true
  optiters=3
[end]

[job]
  type=transient
  debug=dep
[end]