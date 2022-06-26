[mesh]
  type=gmsh
  file=sample.msh
[end]


[dofs]
name=d ux uy
[end]

[elmts]
  [myfracture]
    type=miehefrac
    dofs=d ux uy
    mate=myfracmate
  [end]
[end]

[mates]
  [myfracmate]
    type=miehefracmate
    params=121.15 80.77 2.7e-3 0.02 1.0e-7
    //     lambda mu    Gc     L     viscosity
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=5.0e-10
  r_abs_tol=2.5e-7
  solver=superlu
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
  interval=5
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-5
  adaptive=true
  optiters=3
  growthfactor=1.1
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=1.0e-4
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom left right top
  [end]
  [load]
    type=dirichlet
    dofs=ux
    value=1.0*t
    boundary=top
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]

[postprocess]
  interval=2
  [area]
    type=area
    side=bottom
  [end]
  [ux]
    type=sideintegral
    dof=ux
    side=top
  [end]
  [sigma_xx]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=1
    jindex=1
    side=top
  [end]
  [sigma_xy]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=1
    jindex=2
    side=top
  [end]
  [sigma_yy]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=2
    jindex=2
    side=top
  [end]
  [strain_xx]
    type=rank2matesideintegral
    rank2mate=strain
    iindex=1
    jindex=1
    side=top
  [end]
  [strain_xy]
    type=rank2matesideintegral
    rank2mate=strain
    iindex=1
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

[job]
  type=transient
  debug=dep
[end]
