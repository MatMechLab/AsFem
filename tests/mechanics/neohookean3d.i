*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=3
  zmax=10.0
  nx=4
  ny=4
  nz=100
  meshtype=hex8
[end]

[dofs]
name=ux uy uz
[end]

[projection]
name=von hydro sig_xx sig_yy sig_zz sig_zy sig_zx sig_xy
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux uy uz
    mate=neohookean
    domain=alldomain
  [end]
[end]

[mates]
  [neohookean]
    type=neohookean
    params=1.0e2           0.3           1           
    //     Youngs modulus  Poisson ratio compressive
  [end]
[end]



[bcs]
  [FixUx]
    type=dirichlet
    dof=ux
    boundary=left
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [FixUz]
    type=dirichlet
    dof=uz
    boundary=back
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=ux
    value=1.0e-2*t
    boundary=right
  [end]
  [loadUz]
    type=dirichlet
    dof=uz
    value=1.0e-1*t
    boundary=front
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  endtime=1.0
  adaptive=true
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
  projection=true
  interval=5
[end]
