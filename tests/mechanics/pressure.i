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

[nonlinearsolver]
  type=newtonls
  stol=1.0e-32
[end]


[bcs]
  [FixUx]
    type=dirichlet
    dof=ux
    boundary=right
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [P2Ux]
    type=pressure
    dof=ux
    value=-0.5
    boundary=surface
  [end]
  [P2Uy]
    type=pressure
    dof=uy
    value=-0.5
    boundary=surface
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]
