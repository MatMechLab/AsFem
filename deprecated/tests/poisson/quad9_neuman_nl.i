*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=2.0
  ymin=0.0
  ymax=2.0
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
    mate=nonlinear
    block=alldomain
	[end]
[end]

[mates]
  [nonlinear]
    type=nonlinearpoisson
  [end]
[end]

[linearsolver]
  type=lu
[end]

[nonlinearsolver]
  type=nr
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
    value=-10.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]