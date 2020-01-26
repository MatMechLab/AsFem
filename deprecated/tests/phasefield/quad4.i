*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  nx=100
  ny=100
  meshtype=quad4
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
    params=1.0 2.5 0.04
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
  debug=dep
  dt=1.0e-6
  endtime=2.0e0
  interval=10
[end]