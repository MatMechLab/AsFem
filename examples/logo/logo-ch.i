[mesh]
  type=gmsh
  file=logo.msh
[end]

[dofs]
name=c mu
[end]

[elmts]
  [mych1]
    type=cahnhilliard
    dofs=c mu
    mate=myf1
    domain=matrix1
  [end]
  [mych2]
    type=cahnhilliard
    dofs=c mu
    mate=myf2
    domain=matrix2
  [end]
  [mych3]
    type=cahnhilliard
    dofs=c mu
    mate=myf3
    domain=inclusion
  [end]
[end]

[mates]
  [myf1]
    type=idealsolution
    params=5.0 1.0 0.4
    //     D   Chi Kappa
  [end]
  [myf2]
    type=idealsolution
    params=5.0 1.0 0.4
    //     D   Chi Kappa
  [end]
  [myf3]
    type=idealsolution
    params=1.0 3.5 0.4
    //     D   Chi Kappa
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
  //solver=superlu
[end]

[ics]
  [ic1]
    type=random
    dof=c
    params=0.6 0.66
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]