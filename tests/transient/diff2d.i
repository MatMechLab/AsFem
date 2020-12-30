// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  nx=50
  ny=50
  meshtype=quad9
[end]

[dofs]
name=c
[end]

[qpoint]
  // for quad9 mesh, the order must>=4 !!!
  type=gauss
  order=4
[end]

[elmts]
  [elmt1]
    type=diffusion
    dofs=c
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=constdiffusion
    params=1.0
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=1.0e-2
[end]

[projection]
name=projc dcx dcy
[end]

[ics]
  [randc]
    type=random
    dof=c
    params=1.0 2.0
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]