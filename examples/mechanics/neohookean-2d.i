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
    type=neohookean
    params=210.0 0.3
    //     E     nu
  [end]
[end]

[bcs]
  [fixux]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=left
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom
  [end]
  [loadinguy]
    type=dirichlet
    dofs=uy
    value=1.0*t
    boundary=top
  [end]
[end]

[nonlinearsolver]
  type=newton
  maxiters=20
  r_rel_tol=1.0e-12
  r_abs_tol=6.0e-7
[end]

[timestepping]
  type=be
  dt=2.0e-5
  dtmax=5.0e-1
  //dtmin=5.0e-3
  time=4.0e0
  optiters=12
  adaptive=true
[end]


[output]
  type=vtu
  interval=5
[end]

[projection]
scalarmate=vonMises
rank2mate=stress
[end]

[postprocess]
  [uy]
    type=sideintegral
    dof=uy
    side=top
  [end]
  [sigma_yy]
    type=rank2matesideintegral
    rank2mate=stress
    iindex=2
    jindex=2
    side=top
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]

