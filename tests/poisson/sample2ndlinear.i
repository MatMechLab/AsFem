*** This is an input file for tesing

[mesh]
  type=gmsh
  file=sample2nd.msh
[end]

[dofs]
name=phi
[end]

[qpoint]
  type=gauss
  order=5
[end]

[elmts]
  [poisson]
	  type=poisson
	  dofs=phi
	  mate=linear
    domain=alldomain
  [end]
[end]

[mates]
  [linear]
    type=constpoisson
    params=0.25e0 0.5
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=left right1 right2
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
