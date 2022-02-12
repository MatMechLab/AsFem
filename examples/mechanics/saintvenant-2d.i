[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad9
[end]

[qpoint]
  type=gauss
  order=4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolids]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=saintvenant
    params=210.0 0.3
    //     E     nu
  [end]
[end]

[bcs]
  [fix]
    type=dirichlet
    dofs=ux uy
    value=0.0
    boundary=bottom
  [end]
  [loading]
    type=dirichlet
    dofs=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[nonlinearsolver]
  type=newtontr
  maxiters=25
  r_rel_tol=1.0e-12
  r_abs_tol=4.0e-7
  //solver=superlu
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=5.0e-2
  time=4.0e-1
  optiters=4
  adaptive=true
[end]


[output]
  type=vtu
  interval=20
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[job]
  type=transient
  debug=dep
[end]

