*** 2d neohookean material

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=4.0
  nx=100
  ny=40
  meshtype=quad4
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
  //[loadux]
  //  type=dirichlet
  //  dof=uy
  //  value=0.1*t
  //  boundary=right
  //[end]
  [shearux]
    type=dirichlet
    dof=ux
    value=0.15*t
    boundary=right
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-3
  dtmax=1.0
  endtime=10.0
  adaptive=true
  projection=true
[end]