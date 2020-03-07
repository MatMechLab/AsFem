*** This is an input file for the compressive neohookean model

[mesh]
  type=asfem
  dim=3
  zmax=10.0
  nx=5
  ny=5
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
    params=210.0 0.3
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
  [loadUz]
    type=dirichlet
    dof=uz
    value=0.1
    boundary=front
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]
