[mesh]
  type=gmsh
  file=asymmnotch.msh
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
    params=12.0   8.0 1.0e-3 0.02 1.0e-6
    //     lambda mu  Gc     L     viscosity
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=80
  r_rel_tol=5.0e-10
  r_abs_tol=2.5e-7
  solver=mumps
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
  dt=1.0e-5
  time=5.2e-1
  adaptive=true
  optiters=3
  growthfactor=1.2
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=1.0e-4
[end]

[bcs]
  [fixux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=leftpoint
  [end]
  [fixuy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=leftpoint rightpoint
  [end]
  [load]
    type=dirichlet
    dof=uy
    value=-1.0*t
    boundary=toppoint
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[postprocess]
  [fx]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=1
    jindex=1
    side=top
  [end]
  [fy]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=2
    jindex=2
    side=top
  [end]
  [ux]
    type=sideintegral
    dof=ux
    side=top
  [end]
  [uy]
    type=sideintegral
    dof=uy
    side=top
  [end]
[end]


[job]
  type=transient
  debug=dep
[end]
