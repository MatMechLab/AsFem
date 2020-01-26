*** This is an input file for tesing

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=100
  meshtype=edge2
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

[projection]
name=phi1 dphi
[end]



[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=true
  projection=true
[end]