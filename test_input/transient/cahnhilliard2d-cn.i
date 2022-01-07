// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=80
  ny=80
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
    type=idealsolution
    params=1.0 2.5 0.005
  [end]
[end]

[timestepping]
  type=cn
  dt=1.0e-5
  time=2.0e-5
  optiters=3
  growthfactor=1.2
  adaptive=false
  dtmax=1.0e1
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
  r_rel_tol=1.0e-12
  r_abs_tol=2.0e-7
  solver=mumps
[end]

[projection]
vectormate=gradc
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