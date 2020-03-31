*** This is an input file for nonlinear poisson equation in 2d

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=100
  ny=100
  meshtype=quad4
[end]

//[projection]
//name=gradphi_x gradphi_y
//[end]

[dofs]
name=Phi
[end]

[elmts]
  [poisson]
	  type=user1
	  dofs=Phi
	  mate=nonlinear
    domain=alldomain
  [end]
[end]

[mates]
  [nonlinear]
    type=user1
  [end]
[end]



[bcs]
  [fixphi]
    type=dirichlet
    dof=Phi
    boundary=left
    value=1.0
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]
