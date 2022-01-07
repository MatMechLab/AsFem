[mesh]
  type=gmsh
  file=rect-yz.msh
  savemesh=true
[end]

[dofs]
name=phi
[end]

[qpoint]
type=gauss
order=3
[end]

[elmts]
  [poisson]
    type=poisson
    dofs=phi
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constpoisson
    params=1.0 1.0
  [end]
[end]

[bcs]
  [fixbottom]
    type=dirichlet
    dofs=phi
    value=0.0
    boundary=bottom
  [end]
  [fixtop]
    type=dirichlet
    dofs=phi
    value=1.0
    boundary=top
  [end]
[end]

[projection]
vectormate=gradu
[end]

[job]
type=static
debug=dep
[end]

