[mesh]
  type=asfem
  dim=3
  xmax=1.0
  ymax=1.0
  zmax=5.0
  nx=5
  ny=5
  nz=25
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
    type=saintvenant
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
  type=newton
  //type=newtonsecant
  //type=newtontr
  //type=ncg
  //type=ngmres
  //type=bfgs
  //type=richardson
  //type=nms
  //type=fas
  maxiters=200
  r_rel_tol=1.0e-12
  r_abs_tol=6.0e-7
  //solver=superlu
[end]

[timestepping]
  type=be
  dt=1.0e-3
  dtmax=5.0e-2
  //dtmin=2.0e-3
  time=4.0e-3
  optiters=7
  adaptive=false
[end]


[output]
  type=vtu
  interval=2
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
rank4mate=jacobian
[end]

[job]
  type=transient
  debug=dep
[end]

