[mesh]
  type=asfem
  dim=2
  xmin=-25
  xmax= 25
  ymin=-25
  ymax= 25
  nx=200
  ny=200
  meshtype=quad9
  //meshtype=quad4
[end]

[qpoint]
  type=gauss
  order=4
[end]

[dofs]
name=v u
[end]

[elmts]
  [mymodel]
    type=wave
    dofs=v u
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=wavemate
    params=6.0e-1   1        2.5
    //     speed   choice   radius
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=25
  r_rel_tol=1.0e-9
  r_abs_tol=5.5e-7
  //solver=cg
[end]

[ics]
  [constd]
    type=const
    dof=v
    params=0.0
  [end]
[end]

[output]
  type=vtu
  interval=2
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e2
  adaptive=true
  optiters=4
  growthfactor=1.1
  cutfactor=0.85
  dtmax=1.0e-1
[end]

[projection]
vectormate=gradu
[end]


//[bcs]
//  [fixux]
//    type=dirichlet
//    dofs=u
//    value=0.0
//    boundary=left
//  [end]
//[end]

[job]
  type=transient
  debug=dep
[end]
