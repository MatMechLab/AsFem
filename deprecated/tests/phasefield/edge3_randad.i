*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=5.0
  nx=1000
  meshtype=edge3
[end]

[dofs]
name=conc mu
[end]

[elmts]
  [ch]
    type=cahnhilliard
	  dofs=conc mu
    mate=FreeEnergy
    block=alldomain
  [end]
[end]

[mates]
  [FreeEnergy]
    type=cahnhilliard
    params=5.0 2.5 0.04
    // params= diffusivity chi kappa(for interface energy)
  [end]
[end]


[ics]
  [rand_c]
    type=random
    dof=conc
    values=0.6 0.66
    block=alldomain
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-5
  endtime=8.0e2
  interval=4
  adaptive=true
[end]