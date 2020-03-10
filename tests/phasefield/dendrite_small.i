*** This is an input file for 2d phase field fracture model

[mesh]
  type=asfem
  dim=2
  xmin= 0.0
  xmax= 2.0
  ymin= 0.0
  ymax= 2.0
  nx=250
  ny=250
  meshtype=quad9
[end]

[dofs]
name=phi T
[end]



[projection]
  name=F dFdphi K dkdtheta T_x T_y
[end]


[elmts]
  [dendrite]
    type=dendrite
    dofs=phi T
    mate=dendrite
    domain=alldomain
  [end]
[end]

[mates]
  [dendrite]
    type=dendrite
    params=500.0 0.01 0.04  4 90.0   1.0      1.8
    //     L     eps  delta J theta0 Conduct  eta
  [end]
[end]

[ics]
  [circle]
    type=circle
    dof=phi
    params=1.0 1.0 0.05 0.05  1.0 0.0
    //     x0  y0  R    rwid vin vout
  [end]
[end]




[timestepping]
  type=be
  dt=5.0e-4
  dtmax=5.0e-2
  dtmin=5.0e-7
  endtime=1.0e2
  adaptive=true
[end]
[nonlinearsolver]
  type=newtonls
  r_abs_tol=5.0e-7
  r_rel_tol=1.5e-8
  maxiters=15
[end]


[job]
  type=transient
  debug=dep
  projection=true
[end]
