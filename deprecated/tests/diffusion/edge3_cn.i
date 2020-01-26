*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=100
  meshtype=edge3
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
    params=1.0
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
  endtime=1.0e-3
  timemethod=cn
[end]