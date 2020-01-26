*** 2d neohookean material

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=4.0
  nx=50
  ny=20
  meshtype=quad9
  // I think we need the reduced order integration or volumetric locking modify
  // otherwise, for nonlinear case, the convergence behavior is very bad!!!
  // in contrast, the linear element works very well.
[end]

[dofs]
name=ux uy
[end]

[projection]
name=von hy sig_xx sig_yy sig_xy str_xx str_yy str_xy
[end]



[elmts]
  [solids]
    type=mechanics
	  dofs=ux uy
    mate=neohkmate
    block=alldomain
  [end]
[end]

[mates]
  [neohkmate]
    type=neohookean
    params=1.0e2 0.25
    // params= E nu(Youngs modulus and poisson ratio)
  [end]
[end]

[nonlinearsolver]
  //type=nrbisectlinesearch
  //type=nrregulalinesearch
  type=nr
  stol=1.0e-4
[end]


[bcs]
  [fix_ux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=left
  [end]
  [fix_uy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=left
  [end]
  [loadux]
    type=dirichlet
    dof=ux
    value=0.1
    boundary=top
  [end]
  [shearux]
    type=dirichlet
    dof=uy
    value=0.15
    boundary=right
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]