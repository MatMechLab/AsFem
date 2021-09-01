[mesh]
  type=asfem
  dim=2
  xmax=4.0
  ymax=4.0
  nx=100
  ny=100
  meshtype=quad4
[end]

[dofs]
name=c mu
[end]

[elmts]
  [mych]
    type=cahnhilliard
    dofs=c mu
    mate=myf
  [end]
[end]

[mates]
  [myf]
    type=doublewellpotential
    params=1.0 0.2 0.8 2.0 0.002
    //     D   ca   cb   h kappa
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=3.0e3
  adaptive=true
  optiters=3
  growthfactor=1.1
  cutfactor=0.85
[end]

[nonlinearsolver]
  type=nr
  solver=superlu
[end]

[ics]
  [ic1]
    type=random
    dof=c
    params=0.5 0.55
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]