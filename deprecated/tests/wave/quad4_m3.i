*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1000.0
  ymin=0.0
  ymax=1000.0
  nx=100
  ny=100
  meshtype=quad4
[end]

[dofs]
name=u v
[end]

[elmts]
  [phi]
    type=wave
	  dofs=u v
    mate=mate3
    block=alldomain
  [end]
[end]

[mates]
  [mate3]
    type=generalwave
    params=3.0 2.5 500.0 500.0 20.0 2000.0
    // first is mode choice, second one is the wave speed
  [end]
[end]

//[ics]
//  [randU]
//    type=random
//    dof=u
//    block=alldomain
//  [end]
//[end]


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
  endtime=5.0e4
  adaptive=true
  solver=bicg
[end]