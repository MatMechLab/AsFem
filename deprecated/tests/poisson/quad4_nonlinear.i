*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  nx=400
  ny=400
  meshtype=quad4
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
    params=1.0
  [end]
[end]

[projection]
name=phi0 dphi_x dphi_y
[end]

[linearsolver]
  type=cg
  maxiters=20000
  tol=1.0e-10
[end]

[nonlinearsolver]
  type=nrlinesearch
  s_tol=1.0e-8
[end]

[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=0.0
  [end]
  [right]
    type=dirichlet
    dof=phi
    boundary=right
    value=0.0
  [end]
  [top]
    type=dirichlet
    dof=phi
    boundary=top
    value=0.0
  [end]
  [bottom]
    type=dirichlet
    dof=phi
    boundary=bottom
    value=0.0
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]