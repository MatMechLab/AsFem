//*** this is the input file for step-2 

[mesh]
  type=asfem
  dim=3
  nx=5
  ny=5
  nz=5
  meshtype=hex8
[end]

[qpoint]
type=gauss
order=4
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

[nonlinearsolver]
  type=newton
  //type=ncg
  maxiters=25
  r_rel_tol=1.0e-18
  r_abs_tol=1.0e-15
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
  debug=dep
[end]