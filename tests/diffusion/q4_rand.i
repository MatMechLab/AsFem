*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin= 0.0
  xmax=10.0
  ymin= 0.0
  ymax=10.0
  nx=100
  ny=100
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
    params=1.5e0
  [end]
[end]

[ics]
  [rand]
    type=random
    dof=c
    params=0.1 0.9
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-6
  dtmax=1.0e-1
  opts=4
  endtime=1.0e2
  adaptive=true
[end]

[job]
  type=transient
  debug=dep
  interval=2
[end]
