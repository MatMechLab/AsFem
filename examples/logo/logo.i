//*** this is the input file for step-2 

[mesh]
  type=gmsh
  file=logo.msh
[end]

[dofs]
name=phi
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=phi
    mate=mymate1
    domain=matrix1
  [end]
  [elmt2]
    type=poisson
    dofs=phi
    mate=mymate2
    domain=matrix2
  [end]
  [elmt3]
    type=poisson
    dofs=phi
    mate=mymate3
    domain=inclusion
  [end]
[end]

[mates]
  [mymate1]
    type=constpoisson
    params=2.0e1 1.0e0
  [end]
  [mymate2]
    type=constpoisson
    params=5.0e-2 2.0e0
  [end]
  [mymate3]
    type=constpoisson
    params=5.0e-2 2.0e0
  [end]
[end]

[bcs]
  [fixleft]
    type=dirichlet
    dof=phi
    value=-0.1e4
    boundary=bottom
  [end]
  [fixright]
    type=dirichlet
    dof=phi
    value=-0.1e4
    boundary=top
  [end]
[end]


[job]
  type=static
[end]