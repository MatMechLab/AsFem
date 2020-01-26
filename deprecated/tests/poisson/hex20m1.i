*** This is an input file for tesing

[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  zmin=0.0
  zmax=5.0
  nx=20
  ny=20
  nz=20
  meshtype=hex20
  savemesh=true
[end]

[dofs]
name=phi
[end]

[elmts]
  [phi]
    type=poisson
    dofs=phi
    mate=const
    block=alldomain
  [end]
[end]

[mates]
  [const]
    type=constpoisson
    params=1.0 1.0
  [end]
[end]


[linearsolver]
  type=lu
  maxiters=50000
  tol=1.0e-8
[end]

[bcs]
  [fix1]
    type=dirichlet
    dof=phi
    boundary=left
    value=0.0
  [end]
  [fix2]
    type=dirichlet
    dof=phi
    boundary=right
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]