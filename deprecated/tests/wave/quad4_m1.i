*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  nx=60
  ny=60
  meshtype=quad4
[end]

[dofs]
name=u v
[end]

[elmts]
  [phi]
    type=wave
	  dofs=u v
    mate=mate2
    block=alldomain
  [end]
[end]

[mates]
  [mate2]
    type=generalwave
    params=1.0 10.0
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
  [leftu]
    type=dirichlet
    dof=u
    boundary=left
    value=0.0
  [end]
  [rightu]
    type=dirichlet
    dof=u
    boundary=right
    value=0.0
  [end]
  [bottomu]
    type=dirichlet
    dof=u
    boundary=bottom
    value=0.0
  [end]
  [topu]
    type=dirichlet
    dof=u
    boundary=top
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