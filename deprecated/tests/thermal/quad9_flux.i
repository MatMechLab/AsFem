*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=1.0
  nx=50
  ny=5
  meshtype=quad9
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
    boundary=right
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