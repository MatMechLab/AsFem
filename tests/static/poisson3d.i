//*** this is the input file for step-2 

[mesh]
  type=asfem
  dim=3
  nx=2
  ny=2
  nz=2
  meshtype=hex8
[end]

[dofs]
name=phi
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=phi
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constpoisson
    params=1.0 1.0e1
  [end]
[end]

[bcs]
  [fixleft]
    type=dirichlet
    dofs=phi
    value=0.1
    boundary=left
  [end]
  [fixright]
    type=dirichlet
    dofs=phi
    value=0.5
    boundary=right
  [end]
[end]


[job]
  type=static
[end]