*** simple diffusion for 1d case

[mesh]
  type=gmsh
  file=rect.msh
[end]

[dofs]
name=ux uy d
[end]

[nonlinearsolver]
  //type=nrregulalinesearch
  //type=nrbisectlinesearch
  type=nr
  stol=1.0e-4
  r_abs_tol=2.5e-7
  e_abs_tol=3.5e-15
  maxiters=500
[end]

[linearsolver]
  //type=umfpack
  type=lu
[end]

[elmts]
  [phi]
    type=linearelasticphasefieldfracture
	  dofs=ux uy d
    mate=myfrac
    block=alldomain
  [end]
[end]

[mates]
  [myfrac]
    type=modifyfrac
    params=100.0 0.3 2.7e-3 0.02 1e-6
    // E,nu-->Young's modulus and poisson ratio
    // Gc---> fracture energy
    // L----> thickness
    // viscosity-->for speed of crack propagation
  [end]
[end]

[ics]
  [constD]
    type=const
    dof=d
    values=1.0
    block=alldomain
  [end]
[end]


[bcs]
  [bottomux]
    type=dirichlet
    dof=ux
    boundary=bottom
    value=0.0
  [end]
  [leftuy]
    type=dirichlet
    dof=uy
    boundary=left
    value=0.0
  [end]
  [rightuy]
    type=dirichlet
    dof=uy
    boundary=right
    value=0.0
  [end]
  [bottomuy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [topuy]
    type=dirichlet
    dof=uy
    boundary=top
    value=0.0
  [end]
  [load]
    type=dirichlet
    dof=ux
    boundary=top
    value=1.0*t
  [end]
[end]

[job]
  type=transient
  debug=dep
  dt=1e-4
  dtmax=1.0e-2
  dtmin=2.0e-5
  endtime=5.0e2
  adaptive=true
  optiters=6
  //timemethod=cn
  //projection=true
[end]