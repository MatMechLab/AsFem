[mesh]
  type=asfem
  dim=3
  zmax=10.0
  nx=5
  ny=5
  nz=50
  meshtype=hex8
[end]


[dofs]
name=c ux uy uz
[end]

[elmts]
  [myfracture]
    type=stressdiffusion
    dofs=c ux uy uz
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=diffneohookeanmate
    params=2.0  0.06   150.0 0.2
    //     D    omega  E     nu
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-10
  r_abs_tol=5.0e-7
[end]

[ics]
  [constd]
    type=const
    dof=c
    params=0.0
  [end]
[end]

[output]
  type=vtu
  interval=1
[end]

[timestepping]
  type=be
  dt=1.0e-4
  time=2.0e2
  adaptive=true
  optiters=4
  growthfactor=1.2
  cutfactor=0.85
  dtmax=5.0e-1
[end]

[projection]
scalarmate=vonMises
[end]


[bcs]
  [fixu]
    type=dirichlet
    dofs=ux uy uz
    value=0.0
    boundary=back
  [end]
  [flux]
    type=neumann
    dofs=c
    value=-0.4
    boundary=front
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
