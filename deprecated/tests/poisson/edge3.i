*** This is an input file for tesing

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=2.0
  nx=100000
  meshtype=edge3
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



[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=1.0
  [end]
  [right]
    type=dirichlet
    dof=phi
    boundary=right
    value=2.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]