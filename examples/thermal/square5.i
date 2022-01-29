[mesh]
  type=gmsh
  file=square5.msh
[end]

[dofs]
name=T
[end]

[elmts]
  [grain1]
    type=user1
    dofs=T
    mate=mymate1
    domain=1
  [end]
  [grain2]
    type=user1
    dofs=T
    mate=mymate2
    domain=2
  [end]
  [grain3]
    type=user1
    dofs=T
    mate=mymate3
    domain=3
  [end]
  [grain4]
    type=user1
    dofs=T
    mate=mymate4
    domain=4
  [end]
  [grain5]
    type=user1
    dofs=T
    mate=mymate5
    domain=5
  [end]
[end]

[mates]
  [mymate1]
    type=user2
    params=1.0 2.0 0.05 1.0
    //     rho Cp  K    Q
  [end]
  [mymate2]
    type=user2
    params=1.0 2.0 0.5 1.0
    //     rho Cp  K    Q
  [end]
  [mymate3]
    type=user2
    params=1.0 2.0 5.5 1.0
    //     rho Cp  K    Q
  [end]
  [mymate4]
    type=user2
    params=1.0 2.0 1.5 0.5
    //     rho Cp  K    Q
  [end]
  [mymate5]
    type=user2
    params=1.0 2.0 10.5 1.5
    //     rho Cp  K    Q
  [end]
[end]

[bcs]
  [flux]
    type=neumann
    dofs=T
    value=-0.01
    boundary=top
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