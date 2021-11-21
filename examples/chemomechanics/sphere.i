[mesh]
  type=gmsh
  file=sphere.msh
[end]


[dofs]
name=c mu ux uy uz
[end]

[elmts]
  [mych]
    type=mechcahnhilliard
    dofs=c mu ux uy uz
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelasticchmate
    params=1.0  1.0     0.002  120.0 0.25 0.075  0.5
    //     D    height  kappa  E     nu   Omega  C0
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
  interval=2
[end]

[timestepping]
  type=be
  dt=5.0e-5
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
    type=const
    dof=c
    params=0.5
  [end]
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=ux
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=uy
  [end]
  [fixuz]
    type=dirichlet
    dofs=uz
    value=0.0
    boundary=uz
  [end]
  [flux]
    type=neumann
    dofs=c
    value=0.01
    boundary=surface
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]