*** This is an input file for tesing

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=2000
  meshtype=edge4
[end]

[dofs]
name=u
[end]

[elmts]
	[phi]
	  type=poisson
	  dofs=u
    block=alldomain
	[end]
[end]



[bcs]
  [left]
    type=dirichlet
    dof=u
    boundary=left
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]