*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=2.0
  nx=200
  ny=40
  meshtype=quad8
[end]

[dofs]
name=conc
[end]

[elmts]
  [phi]
    type=diffusion
	  dofs=conc
    mate=const
    block=alldomain
  [end]
[end]

[mates]
  [const]
    type=constdiffusion
  [end]
[end]



[bcs]
  [right]
    type=neumann
    dof=conc
    boundary=right
    value=-1.0
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-5
  endtime=1.0e1
  interval=2
  adaptive=true
[end]