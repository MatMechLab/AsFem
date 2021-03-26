[mesh]
  type=gmsh
  file=bending.msh
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
    params=12.0   8.0 5.0e-4 0.03 1.0e-6
    //     lambda mu  Gc     L     viscosity
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
  time=2.0e-1
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
[end]


[job]
  type=transient
  debug=dep
[end]
