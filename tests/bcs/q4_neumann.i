*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmax= 1.0
  ymax= 1.0
  nx=10
  ny=10
  meshtype=quad4
[end]

[dofs]
	name=c
[end]


[elmts]
  [poisson]
	  type=diffusion
	  dofs=c
	  mate=linear
    domain=alldomain
  [end]
[end]

[mates]
  [linear]
    type=constdiffusion
    params=1.0e0
  [end]
[end]

[bcs]
  [mybc1]
    type=dirichlet
    dofs=c
    boundary=left
    value=0.0
  [end]
  [mybc2]
    type=neumann
    dofs=c
    boundary=right
    value=0.01
  [end]
[end]

[ics]
  [rand]
    type=random
    dof=c
    params=0.0 1.0
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=2.0e-5
  optiters=4
  time=2.0e-5
  adaptive=false
[end]

[job]
  type=transient
  debug=dep
[end]
