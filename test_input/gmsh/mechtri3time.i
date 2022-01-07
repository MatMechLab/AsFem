[mesh]
  type=gmsh
  file=sample.msh
  savemesh=true
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-12
  r_abs_tol=2.6e-7
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
    dofs=ux
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom left right top
  [end]
  [load]
    type=dirichlet
    dofs=ux
    value=t*0.1
    boundary=top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
