[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolids]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
    //     E     nu
  [end]
[end]

[bcs]
  [fix]
    type=dirichlet
    dofs=ux uy
    value=0.0
    boundary=bottom
  [end]
  [loading]
    type=cyclicdirichlet
    dofs=uy
    params= 0.0 0.0 1.0 1.0
    boundary=top
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=1.0e-2
  time=2.5e-1
  optiters=3
  adaptive=true
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[job]
  type=transient
  debug=dep
[end]

