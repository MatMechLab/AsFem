[mesh]
  type=asfem
  dim=3
  xmax=0.1
  ymax=0.1
  zmax=1.0
  nx=5
  ny=5
  nz=50
  meshtype=hex8
[end]

[qpoint]
  type=gauss
  order=4
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
[end]

[nonlinearsolver]
  type=newtonsecant
  maxiters=20
  r_rel_tol=1.0e-12
  r_abs_tol=6.0e-7
  //solver=superlu
[end]

[output]
  type=vtu
  interval=20
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=1.0e-2
  time=5.0e-1
  optiters=5
  adaptive=true
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[job]
  type=transient
  debug=dep
[end]
