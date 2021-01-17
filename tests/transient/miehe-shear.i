[mesh]
  type=gmsh
  file=sample.msh
  savemesh=true
[end]

[dofs]
name=d ux uy
[end]

[elmts]
  [myfracture]
    type=miehefrac
    dofs=d ux uy
    mate=myfracmate
    domain=block
  [end]
[end]

[mates]
  [myfracmate]
    type=miehefracmate
    params=121.15 80.77 2.7e-3 0.015 1.0e-6
    //     lambda mu    Gc     L     viscosity
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-9
  r_abs_tol=4.6e-7
  //solver=superlu
[end]

[timestepping]
  type=be
  dt=1.0e-6
  time=2.5e-2
  adaptive=true
  optiters=3
  growthfactor=1.1
  cutfactor=0.85
  dtmin=1.0e-10
  dtmax=1.0e1
[end]

[bcs]
  [fixux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom left right top
  [end]
  [load]
    type=dirichlet
    dof=ux
    value=1.0*t
    boundary=top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
