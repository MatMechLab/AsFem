[mesh]
  type=asfem
  dim=3
  xmax=1.0
  ymax=1.0
  zmax=10.0
  nx=5
  ny=5
  nz=50
  meshtype=hex8
[end]

[dofs]
name=c
[end]

[elmts]
  [mydiffusion]
    type=diffusion
    dofs=c
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constdiffusion
    params=1.0e1
  [end]
[end]

[bcs]
  [flux]
    type=neumann
    dofs=c
    value=-0.5
    boundary=front
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-3
  time=1.0e0
[end]

[job]
  type=transient
  debug=dep
[end]