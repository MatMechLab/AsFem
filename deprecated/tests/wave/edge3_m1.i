*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=10.0
  nx=2000
  meshtype=edge3
[end]

[dofs]
name=u v
[end]

[elmts]
  [phi]
    type=wave
	  dofs=u v
    mate=mate1
    block=alldomain
  [end]
[end]

[mates]
  [mate1]
    type=generalwave
    params=1.0 1.0
    // first is mode choice, second one is the wave speed
  [end]
[end]

[ics]
  [randU]
    type=random
    dof=u
    block=alldomain
  [end]
[end]


[bcs]
  [left]
    type=dirichlet
    dof=v
    boundary=left
    value=0.0
  [end]
  [right]
    type=dirichlet
    dof=v
    boundary=right
    value=0.0
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=2.0e-6
  dtmax=5.0
  endtime=1.0e2
  adaptive=true
[end]