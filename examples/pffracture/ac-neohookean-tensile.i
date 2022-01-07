[mesh]
  type=gmsh
  file=sample.msh
[end]


[dofs]
name=d ux uy
[end]

[elmts]
  [myfracture]
    type=allencahnfrac
    dofs=d ux uy
    mate=myfracmate
  [end]
[end]

[mates]
  [myfracmate]
    type=neohookeanpffracturemate
    params=121.15 80.77 2.7e-3 0.012 1.0e-6
    //     lambda mu    Gc     L     viscosity(M=1.0/viscosity)
  [end]
[end]

[nonlinearsolver]
  type=newton
  maxiters=25
  r_rel_tol=1.0e-10
  r_abs_tol=5.5e-7
  //solver=superlu
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
  interval=10
[end]

[timestepping]
  type=be
  dt=1.0e-6
  time=2.0e-1
  adaptive=true
  optiters=7
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
    dof=uy
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
    dofs=ux uy
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=top
  [end]
  [load]
    type=dirichlet
    dofs=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
