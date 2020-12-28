// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  nx=100
  ny=100
  meshtype=quad9
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
    dof=u
    value=0.1
    boundary=left
  [end]
  [fixright]
    type=dirichlet
    dof=u
    value=0.5
     boundary=right
  [end]
[end]