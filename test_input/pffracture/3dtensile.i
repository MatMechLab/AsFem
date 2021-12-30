[mesh]
  type=gmsh
  file=sample3d.msh
[end]


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
    type=miehefracmate
    params=12.0   8.0 5.0e-4 0.2 1.0e-6
    //     lambda mu  Gc     L   viscosity
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
  //interval=5
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-5
  adaptive=true
  optiters=3
  growthfactor=1.2
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=2.5e-4
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
    boundary=bottom
  [end]
  [fixuz]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom
  [end]
  [load]
    type=dirichlet
    dofs=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[projection]
scalarmate=vonMises
[end]

[postprocess]
  [ux]
    type=sideintegral
    dof=ux
    side=top
  [end]
  [uy]
    type=sideintegral
    dof=uy
    side=top topload
  [end]
[end]


[job]
  type=transient
  debug=dep
[end]
