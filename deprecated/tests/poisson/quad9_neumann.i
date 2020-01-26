*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  nx=100
  ny=100
  meshtype=quad9
[end]

[dofs]
name=phi
[end]

[elmts]
	[phi]
	  type=poisson
	  dofs=phi
    mate=const
    block=alldomain
	[end]
[end]

[mates]
  [const]
    type=constpoisson
    params=1.0 2.0
  [end]
[end]

[linearsolver]
  type=cg
[end]

[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=0.0
  [end]
  [top]
    type=neumann
    dof=phi
    boundary=right
    value=-2.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]