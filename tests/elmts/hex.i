// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=3
  nx=10
  ny=10
  nz=10
  meshtype=hex27
  printmesh=true
  savemesh=true
[end]

[dofs]
name=disp_x disp_y
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=disp_x disp_y
  [end]
[end]