[mesh]
  type=gmsh
  file=sample.msh
[end]


[dofs]
name=c ux uy
[end]

[elmts]
  [myelmt]
    type=user1
    dofs=c ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=user1
    params=210.0 0.3 1.0
    //     E     nu  D
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-9
  r_abs_tol=4.6e-7
  //solver=superlu
[end]

[ics]
  [constd]
    type=const
    dof=c
    params=0.0
  [end]
  [constux]
    type=const
    dof=ux
    params=0.1
  [end]
  [constuy]
    type=const
    dof=uy
    params=0.05
  [end]
[end]


[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-5
  adaptive=true
  optiters=3
  growthfactor=1.1
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=1.0e1
[end]



[bcs]
  [fixux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom left right top
  [end]
  [load]
    type=dirichlet
    dof=ux
    value=0.01*t
    boundary=top
  [end]
  [flux1]
    type=dirichlet
    dof=c
    value=0.0
    boundary=top
  [end]
  [flux2]
    type=dirichlet
    dof=c
    value=0.1
    boundary=bottom
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
