*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=10.0
  nx=5
  ny=5
  nz=100
  meshtype=hex8
[end]

[dofs]
name=T
[end]

[elmts]
  [phi]
    type=diffusion
	  dofs=T
    mate=const
    block=alldomain
  [end]
[end]

[mates]
  [const]
    type=constdiffusion
    params=1.0
  [end]
[end]

[ics]
  [const]
    type=const
    dof=T
    values=100.0
    block=alldomain
  [end]
[end]


[bcs]
  [right]
    type=neumann
    dof=T
    boundary=front
    value=-5.0
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=5.0e-3
  endtime=1.0
  solver=bicg
  //adaptive=true
[end]