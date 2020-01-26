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
    mate=user1
    block=alldomain
  [end]
[end]

[mates]
  [user1]
    type=user1
    params=2.0
  [end]
[end]

[ics]
  [rand]
    type=random
    dof=conc
    values=0.25 0.35
    block=alldomain
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
  endtime=1.0e-2
  interval=10
[end]