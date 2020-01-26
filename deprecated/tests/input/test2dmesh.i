*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=2.0
  ymin=0.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=u
[end]

[elmts]
	[phi]
	  type=poisson1d
	  dofs=u
    block=alldomain
	[end]
[end]

[bcs]
  [left]
    type=dirichlet
    dof=u
    boundary=left
  [end]
[end]