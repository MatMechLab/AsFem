*** This is an input file for tesing

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=2500
  meshtype=edge4
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
    params=1.0 2.0
  [end]
[end]

[linearsolver]
  type=cg
  maxiters=20000
  tol=1.0e-10
  restart=1500
[end]

[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=1.0
  [end]
  [right]
    type=neumann
    dof=phi
    boundary=right
    value=2.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]