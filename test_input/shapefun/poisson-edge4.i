[mesh]
  type=asfem
  dim=1
  nx=10
  meshtype=edge4
[end]

[dofs]
name=phi
[end]

[elmts]
  [myp]
    type=poisson
    dofs=phi
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constpoisson
    params=0.01 1000
  [end]
[end]

[bcs]
  [left]
    type=dirichlet
    dofs=phi
    value=0.0
    boundary=left
  [end]
  [right]
    type=dirichlet
    dofs=phi
    value=1.0
    boundary=right
  [end]
[end]


[projection]
vectormate=gradu
[end]


[job]
type=static
debug=dep
[end]