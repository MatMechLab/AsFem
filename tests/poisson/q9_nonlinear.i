*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=-5.0
  xmax= 5.0
  ymin=-5.0
  ymax= 5.0
  nx=100
  ny=100
  meshtype=quad9
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
    params=1.5e0 0.5
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=left right bottom top
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
