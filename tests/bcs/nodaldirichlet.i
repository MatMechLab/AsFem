// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  xmax=10.0
  ymax=2.0
  nx=200
  ny=40
  meshtype=quad4
[end]

[dofs]
name=c
[end]

[elmts]
  [elmt1]
    //type=helper
    type=diffusion
    dofs=c
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=constdiffusion
    params=1.0e2
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-5
  adaptive=false
  optiters=3
[end]

[projection]
vectormate=gradc
[end]


[bcs]
  [fix]
    //type=helper
    type=nodaldirichlet
    dof=c
    value=1.0
    boundary=left bottom
  [end]
  [fix1]
    type=nodaldirichlet
    dof=c
    value=2.0
    boundary=left bottom
  [end]
[end]

[job]
  //type=helper
  type=transient
  debug=dep
[end]
