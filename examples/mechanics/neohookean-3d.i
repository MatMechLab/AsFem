[mesh]
  type=asfem
  dim=3
  xmax=1.0
  ymax=1.0
  zmax=4.0
  nx=10
  ny=10
  nz=50
  meshtype=hex8
[end]

[qpoint]
  type=gauss
  order=2
[end]

[dofs]
name=ux uy uz
[end]

[elmts]
  [mysolids]
    type=mechanics
    dofs=ux uy uz
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=neohookean
    params=210.0 0.3
    //     E     nu
  [end]
[end]

[bcs]
  [fix]
    type=dirichlet
    dofs=ux uy uz
    value=0.0
    boundary=back
  [end]
  [loading]
    type=dirichlet
    dofs=uz
    value=1.0*t
    boundary=front
  [end]
  [loadinguy]
    type=dirichlet
    dofs=uy
    value=0.5*t
    boundary=front
  [end]
[end]

[nonlinearsolver]
  type=newtonsecant
  maxiters=200
  r_rel_tol=1.0e-12
  r_abs_tol=6.0e-7
  solver=superlu
[end]

[timestepping]
  type=be
  dt=2.0e-2
  dtmax=5.0e-2
  //dtmin=5.0e-3
  time=4.0e0
  optiters=7
  adaptive=false
[end]


[output]
  type=vtu
  interval=1
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[job]
  type=transient
  debug=dep
[end]

