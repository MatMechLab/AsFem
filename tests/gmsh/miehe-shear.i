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
    //type=helper
    type=const
    dof=d
    params=0.0
  [end]
[end]

[output]
  //type=helper
  type=vtu
  interval=3
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

[projection]
name=reacforce_x reacforce_y
scalarmate=vonMises
rank2mate=stress
[end]

[postprocess]
  [area]
    type=area
    side=bottom
  [end]
  [fx]
    type=projvariablesideintegral
    projvariable=reacforce_x
    side=top
  [end]
  [fy]
    type=projvariablesideintegral
    projvariable=reacforce_y
    side=top
  [end]
  [ux]
    type=sideintegral
    dof=ux
    side=top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
