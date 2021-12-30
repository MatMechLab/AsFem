// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=3
  nx=10
  ny=10
  nz=10
  meshtype=hex27
  savemesh=true
[end]

[dofs]
name=u
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=u
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=constpoisson
    params=1.0 1.0e1
  [end]
[end]

[bcs]
  [fixleft]
    type=dirichlet
    dofs=u
    value=0.1
    boundary=left
  [end]
  [fixright]
    type=dirichlet
    dofs=u
    value=0.5
     boundary=right
  [end]
[end]

[qpoint]
  type=gauss
  order=4
  bcorder=2
[end]

[projection]
vectormate=gradu
[end]

[job]
  type=static
  debug=dep
[end]