[mesh]
  type=gmsh
  file=holes3d.msh
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
    type=neohookeanpffracturemate
    params=121.15 80.77 2.7e-3 0.015 1.0e-6
    //     lambda mu    Gc     L     viscosity
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-10
  r_abs_tol=4.5e-7
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
  interval=20
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
[end]



[bcs]
  [fixux]
    type=dirichlet
    dofs=ux uy uz
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
