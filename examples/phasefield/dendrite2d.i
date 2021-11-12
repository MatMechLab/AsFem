[mesh]
  type=asfem
  dim=2
  xmax=10.0
  ymax=10.0
  nx=500
  ny=500
  meshtype=quad4
[end]

[qpoint]
  type=gauss
  order=4
[end]


[dofs]
name=eta T
[end]

[elmts]
  [mydendrite]
    type=kobayashi
    dofs=eta T
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=kobayashimate
    params=3.0e3  0.02 0.1   4.0  -1.8
    //     L      K0   delta N    latent-heat
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-10
  r_abs_tol=4.5e-7
  //solver=superlu
[end]

[ics]
  [const]
    type=const
    dof=T
    params=1.0
  [end]
  [nuclei]
    type=smoothcircle
    dof=eta
    params=0.0 5.0 0.1 0.1  1.0
    //     x0  y0  r    dr    value
    domain=alldomain
  [end]
[end]

[output]
  type=vtu
  interval=1
[end]

[timestepping]
  type=be
  dt=1.0e-4
  time=2.0e-1
  adaptive=true
  optiters=6
  growthfactor=1.2
  cutfactor=0.85
  dtmax=1.0e-1
[end]

[projection]
vectormate=GradEta
[end]


[job]
  type=transient
  debug=dep
[end]
