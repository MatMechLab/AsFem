*** This is an input file for tesing

[mesh]
  type=abaqus
  file=quad8_fine.inp
[end]

[dofs]
name=phi
[end]

[elmts]
  [poisson]
	  type=poisson
	  dofs=phi
	  mate=const
    domain=alldomain
  [end]
[end]

[mates]
  [const]
    type=constpoisson
    params=1.5e0 0.5
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=bottom top
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
