*** This is an input file for tesing

[mesh]
  type=gmsh
  file=sample.msh
[end]

[dofs]
name=phi
[end]

[elmts]
  [poisson]
	  type=poisson
	  dofs=phi
	  mate=nonlinear
    domain=alldomain
  [end]
[end]

[mates]
  [nonlinear]
    type=nlpoisson
    params=0.25e0 0.5
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=left bottom1 bottom2 right1 right2 top
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
