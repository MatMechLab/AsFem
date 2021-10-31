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
    type=stressdecompositionmate
    params=121.15 80.77 2.7e-3 0.012  1.0e-6     0.0 0.0 0.0 1  
    //     lambda mu    Gc     L      viscosity  t1  t2  t3  UseHist
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
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
  interval=10
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
  dtmax=1.0e-3
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux uy
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
rank2mate=stress
[end]


[job]
  type=transient
  debug=dep
[end]
