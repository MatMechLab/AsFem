*** 3d neohookean material

[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  zmin=0.0
  zmax=60.0
  nx=5
  ny=5
  nz=80
  meshtype=hex8
[end]

[dofs]
name=disp_x disp_y disp_z
[end]

[projection]
name=von hy sig_xx sig_yy sig_zz sig_zy sig_zx sig_xy
[end]

[elmts]
  [solids]
    type=mechanics
	  dofs=disp_x disp_y disp_z
    mate=neohkmate
    block=alldomain
  [end]
[end]

[mates]
  [neohkmate]
    type=neohookean
    params=1.0e3 0.40
    // params= E nu(Youngs modulus and poisson ratio)
  [end]
[end]


[bcs]
  [fix_ux]
    type=dirichlet
    dof=disp_x
    value=0.0
    boundary=back
  [end]
  [fix_uy]
    type=dirichlet
    dof=disp_y
    value=0.0
    boundary=back
  [end]
  [fix_uz]
    type=dirichlet
    dof=disp_z
    value=0.0
    boundary=back
  [end]
  [loadux]
    type=dirichlet
    dof=disp_y
    value=0.1*t
    boundary=front
  [end]
  [shearux]
    type=dirichlet
    dof=disp_z
    value=-0.2*t
    boundary=front
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-3
  dtmax=1.0e0
  endtime=1.0e3
  optiters=4
  adaptive=true
  projection=true
[end]