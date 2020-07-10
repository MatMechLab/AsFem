// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=1
  nx=100
  meshtype=edge2
  printmesh=true
  savemesh=true
[end]

[dofs]
name=u v
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=u v
  [end]
  [elmt2]
    type=diffusion
    dofs=u v
  [end]
[end]