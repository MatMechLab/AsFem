*** This is an input file for tesing

[mesh]
  type=asfem
  dim=3
  xmax= 1.0
  ymax= 1.0
  zmax=10.0
  nx=5
  ny=5
  nz=50
  meshtype=hex27
[end]

[qpoint]
  type=gauss
  order=4
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
    params=5.0e0
  [end]
[end]

[bcs]
  [mybc2]
    type=neumann
    dofs=c
    boundary=front
    value=-0.01
  [end]
[end]



[timestepping]
  type=be
  dt=1.0e-5
  dtmax=1.0e-1
  optiters=4
  time=2.0e-2
  adaptive=true
[end]

[job]
  type=transient
  debug=dep
[end]
