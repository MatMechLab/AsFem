[mesh]
  type=asfem
  dim=3
  xmax=0.8
  ymax=0.8
  zmax=8
  nx=4
  ny=4
  nz=100
  meshtype=hex8
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
    params=5.0  1.0     0.002  120.0 0.25 0.075  0.12
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
    params=0.12
  [end]
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=back
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=back
  [end]
  [fixuz]
    type=dirichlet
    dofs=uz
    value=0.0
    boundary=back
  [end]
  [flux]
    type=neumann
    dofs=c
    value=-0.01
    boundary=front
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]