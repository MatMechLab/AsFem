// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
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
    params=1.0 2.5 0.002
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-4
  time=1.0e0
  optiters=3
  growthfactor=1.2
  adaptive=true
  dtmin=1.0e-8
  dtmax=1.0e1
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-8
  r_abs_tol=1.0e-7
  solver=mumps
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

[postprocess]
  [totalc]
    type=elementalintegral
    dof=c
    domain=alldomain
  [end]
  [area]
    type=volume
    domain=alldomain
  [end]
  [nodec]
    type=nodevalue
    nodeid=315
  [end]
  [elmtc]
    type=elementvalue
    elmtid=20
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]