[mesh]
  type=gmsh
  file=line3d.msh
  savemesh=true
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
    params=1 1
  [end]
[end]

[bcs]
  [left]
    type=dirichlet
    dof=phi
    value=0.0
    boundary=left
  [end]
  [right]
    type=dirichlet
    dof=phi
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