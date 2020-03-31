*** This is an input file for tesing

[mesh]
  type=gmsh
  file=sample2nd.msh
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
    params=0.25e1 1.5
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=left
    value=0.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
