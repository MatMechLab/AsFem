*** simple diffusion for 1d case

[mesh]
  type=gmsh
  file=circle.msh
  savemesh=true
[end]

[dofs]
name=conc mu
[end]

[projection]
name=F c chemical cx cy
[end]

[elmts]
  [ch]
    type=cahnhilliard
	  dofs=conc mu
    mate=FreeEnergy
  [end]
[end]

[linearsolver]
  type=lu
[end]

[mates]
  [FreeEnergy]
    type=cahnhilliard
    params=1.5 2.5 0.025
    // params= diffusivity chi kappa(for interface energy)
  [end]
[end]


[ics]
  [rand_c]
    type=random
    dof=conc
    values=0.45 0.55
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=2.0e-6
  dtmax=0.05
  endtime=5.0e2
  //interval=2
  adaptive=true
  projection=true
[end]