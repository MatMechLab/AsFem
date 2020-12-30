// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  xmax=8.0
  ymax=8.0
  nx=60
  ny=60
  meshtype=quad9
[end]

[dofs]
name=c mu
[end]

[qpoint]
  // for quad9 mesh, the order must>=4 !!!
  type=gauss
  order=4
[end]

[elmts]
  [elmt1]
    type=cahnhilliard
    dofs=c mu
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=cahnhilliard
    params=1.0 2.5 0.03
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-2
  time=1.0e1
[end]

[projection]
name=projc dcx dcy
[end]

[ics]
  [randc]
    type=random
    dof=c
    params=0.6 0.63
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]