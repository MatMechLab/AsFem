*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=2.0
  ymin=0.0
  ymax=2.0
  zmin=0.0
  zmax=10.0
  nx=5
  ny=5
  nz=120
  meshtype=hex8
[end]

[dofs]
name=ux uy uz c
[end]

[elmts]
  [phi]
    type=thermalmechanics
	  dofs=ux uy uz c
    mate=fullcouple
    block=alldomain
  [end]
[end]

[mates]
  [fullcouple]
    type=thermelastic
    params=1.0 0.1 0.08 100.0 0.25 1.0 2.0 0.5 0.0 0.0 0.0
    // D,Cref,Omega,E,nu, eigen strain
    // D->diffusivity(conductivity)
    // Cref->reference concentration or temperature
    // E,nu-->Young's modulus and poisson ratio
    // eigen strain(in voigt notation, even in 2D, 6-components are requried!!!)
  [end]
[end]

[ics]
  [randU]
    type=const
    dof=c
    values=0.1
    block=alldomain
  [end]
[end]

[linearsolver]
  type=cg
  maxiters=50000
[end]

[bcs]
  [leftux]
    type=dirichlet
    dof=ux
    boundary=left
    value=0.0
  [end]
  [bottomuy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [backuz]
    type=dirichlet
    dof=uz
    boundary=back
    value=0.0
  [end]
  [flux]
    type=neumann
    dof=c
    boundary=front
    value=-0.01
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-6
  dtmax=0.1
  endtime=5.0e1
  adaptive=true
[end]