[mesh]
  type=asfem
  dim=3
  xmax=0.04
  ymax=0.04
  zmax=0.4
  nx=10
  ny=10
  nz=100
  meshtype=hex8
[end]

//[qpoint]
//  type=gauss
//  order=4
//[end]

[dofs]
name=d ux uy uz
[end]

[elmts]
  [myfracture]
    type=miehefrac
    dofs=d ux uy uz
    mate=myfracmate
  [end]
[end]

[mates]
  [myfracmate]
    type=neohookeanpffracturemate
    params=121.15 80.77 2.7e-3 0.002  1.0e-6
    //     lambda mu    Gc     L      viscosity
  [end]
[end]

[nonlinearsolver]
  type=newtonsecant
  maxiters=20
  r_rel_tol=1.0e-10
  r_abs_tol=5.5e-7
  //solver=superlu
  //solver=cg
[end]

[ics]
  [constd]
    type=const
    dof=d
    params=0.0
  [end]
[end]

[output]
  type=vtu
  interval=1
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-1
  adaptive=true
  optiters=5
  growthfactor=1.2
  cutfactor=0.85
  dtmax=1.0e-4
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]

[postprocess]
  [uy]
    type=sideintegral
    dof=uz
    side=top
  [end]
  [sigma_yy]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=2
    jindex=2
    side=top
  [end]
  [strain_yy]
    type=rank2matesideintegral
    rank2mate=strain
    iindex=2
    jindex=2
    side=top
  [end]
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux uy uz
    value=0.0
    boundary=back
  [end]
  [fixuy]
    type=dirichlet
    dofs=ux uy
    value=0.2*t
    boundary=front
  [end]
  [load]
    type=dirichlet
    dofs=uz
    value=1.0*t
    boundary=front
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
