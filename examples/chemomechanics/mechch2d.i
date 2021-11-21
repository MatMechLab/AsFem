[mesh]
  type=asfem
  dim=2
  xmax=2
  ymax=2
  nx=50
  ny=50
  meshtype=quad4
[end]


[dofs]
name=c mu ux uy
[end]

[elmts]
  [mych]
    type=mechcahnhilliard
    dofs=c mu ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelasticchmate
    params=1.0  1.0     0.01  120.0 0.25 0.06  0.5
    //     D    height  kappa E     nu   Omega C0
    // D     : diffusivity
    // height: energy barrier height
    // E     : Youngs modulus
    // nu    : poisson ratio
    // kappa : interface thickness parameter
    // Omega : partial molar volume
    // C0    : reference concentration
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-10
  r_abs_tol=5.0e-7
[end]


[output]
  type=vtu
  interval=5
[end]

[timestepping]
  type=be
  dt=5.0e-4
  time=1.0e2
  adaptive=true
  optiters=4
  growthfactor=1.2
  cutfactor=0.85
  dtmax=2.0e-1
[end]

[projection]
scalarmate=vonMises HyStress
[end]

[ics]
  [constd]
    type=random
    dof=c
    params=0.45 0.55
  [end]
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
[end]

[job]
  type=transient
  debug=dep
[end]