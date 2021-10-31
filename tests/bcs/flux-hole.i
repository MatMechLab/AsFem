*** This is an input file for tesing

[mesh]
  type=gmsh
  file=hole.msh
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
  [mybc2]
    type=neumann
    dofs=c
    boundary=surface
    value=-0.01
  [end]
[end]


[timestepping]
  type=be
  dt=1.0e-5
  dtmax=4.0e-5
  optiters=4
  time=2.0e-2
  adaptive=true
[end]

[job]
  type=transient
  debug=dep
[end]
