[mesh]
  type=asfem
  dim=3
  xmax=2.0
  ymax=2.0
  zmax=10.0
  nx=10
  ny=10
  nz=50
  meshtype=hex8
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
    type=linearelastic
    params=210.0 0.3  
    //     E     nu
  [end]
[end]

[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux uy uz
    boundary=back
    value=0.0
  [end]
  [loadUz]
    type=dirichlet
    dofs=uz
    value=1.0*t
    boundary=front
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
  r_rel_tol=5.0e-8
  r_abs_tol=4.5e-7
  //solver=cg
  //solver=mumps
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=1.0e-2
  dtmax=1.0
  adaptive=true
  optiters=3
[end]

[job]
  type=transient
  debug=dep
[end]