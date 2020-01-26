*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=150
  ny=30
  meshtype=quad4
[end]

[dofs]
name=ux uy c
[end]

[elmts]
  [phi]
    type=thermalmechanics
	  dofs=ux uy c
    mate=fullcouple
    block=alldomain
  [end]
[end]

[mates]
  [fullcouple]
    type=thermelastic
    params=1.0 0.1 0.08 100.0 0.25 0.5 1.0 1.0 0.0 0.0 0.1
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
  [fluxright]
    type=neumann
    dof=c
    boundary=right
    value=-0.01
  [end]
  [fluxtop]
    type=neumann
    dof=c
    boundary=top
    value=-0.01
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-6
  dtmax=1.0
  endtime=5.0e1
  adaptive=true
  solver=bicg
  projection=true
[end]