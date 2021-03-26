[mesh]
  type=gmsh
  file=platewithhole.msh
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
  [fixright]
    type=dirichlet
    dof=phi
    value=2.0
    boundary=right
  [end]
  [fixleft]
    type=dirichlet
    dof=phi
    value=1.0
    boundary=left
  [end]
  [fixpoint]
    type=nodaldirichlet
    dof=phi
    value=-5.0
    boundary=bottompoint
  [end]
[end]

[job]
type=static
debug=dep
[end]

