*** This is an input file for pressure loading test

[mesh]
  type=gmsh
  file=platehole.msh
[end]

[dofs]
name=ux uy
[end]

[projection]
name=von hydro sig_xx sig_yy sig_xy
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux uy
    mate=elastic
    domain=alldomain
  [end]
[end]

[mates]
  [elastic]
    type=linearelastic
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
    boundary=top
    value=0.0
  [end]
  [P2Ux]
    type=pressure
    dof=ux
    value=2.5*t
    boundary=surface
  [end]
  [P2Uy]
    type=pressure
    dof=uy
    value=2.5*t
    boundary=surface
  [end]
[end]


[timestepping]
  type=be
  dt=1.0e-3
  endtime=0.4
  adaptive=true
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
  projection=true
  //interval=5
[end]
