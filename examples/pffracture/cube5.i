[mesh]
  type=gmsh
  file=cube5.msh
  savemesh=true
[end]

[dofs]
name=d ux uy uz
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
  dt=2.0e-5
  time=4.0e-5
  adaptive=true
  optiters=4
  growthfactor=1.2
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=2.0e-3
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=left
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom
  [end]
  [fixuz]
    type=dirichlet
    dofs=uz
    value=0.0
    boundary=back
  [end]
  [load]
    type=dirichlet
    dofs=uz
    value=1.0*t
    boundary=front
  [end]
[end]

[projection]
scalarmate=vonMises
[end]

[job]
  type=transient
  debug=dep
[end]


//***************************************************
[elmts]
  [myelmt1]
    type=miehefrac
    dofs=d ux uy uz
    mate=mymate1
    domain=1
  [end]
  [myelmt2]
    type=miehefrac
    dofs=d ux uy uz
    mate=mymate2
    domain=2
  [end]
  [myelmt3]
    type=miehefrac
    dofs=d ux uy uz
    mate=mymate3
    domain=3
  [end]
  [myelmt4]
    type=miehefrac
    dofs=d ux uy uz
    mate=mymate4
    domain=4
  [end]
  [myelmt5]
    type=miehefrac
    dofs=d ux uy uz
    mate=mymate5
    domain=5
  [end]
[end]


[mates]
  [mymate1]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  5.7787e+01  0.0000e+00  0.0000e+00
  [end]
  [mymate2]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.2437e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate3]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.7698e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate4]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  1.4904e+02  0.0000e+00  0.0000e+00
  [end]
  [mymate5]
    type=stressdecompositionmate
    params= 2.8269e+02  1.4807e+02  1.4807e+02  1.4134e+02  6.0575e+01  1.4134e+02  4.0385e+01  4.0385e+01  8.0770e+01  2.7000e-03  1.2000e-02  1.0000e-06  8.4469e+01  0.0000e+00  0.0000e+00
  [end]
[end]
