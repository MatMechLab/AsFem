*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  nx=250
  ny=250
  meshtype=quad4
[end]

[dofs]
name=phi
[end]

[elmts]
	[phi]
	  type=poisson
	  dofs=phi
    block=alldomain
    mate=const
	[end]
[end]

[mates]
  [const]
    type=constpoisson
    params=1.0 2.0
  [end]
[end]

[projection]
name=phi0 dphi_x dphi_y
[end]

[linearsolver]
  type=cg
  maxiters=2000
  tol=1.0e-8
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
  projection=true
[end]