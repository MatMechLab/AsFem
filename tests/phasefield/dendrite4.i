*** This is an input file for 2d phase field fracture model

[mesh]
  type=asfem
  dim=2
  xmin= 0.0
  xmax= 4.0
  ymin= 0.0
  ymax= 4.0
  nx=300
  ny=300
  meshtype=quad4
[end]

[dofs]
name=phi T
[end]


[qpoint]
  type=gauss
  order=4
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
    params=3333.33 0.03  0.08  4 90.0   1.0      1.8
    //     L       eps   delta J theta0 Conduct  eta
  [end]
[end]

[ics]
  [circle]
    type=circle
    dof=phi
    params=2.0 2.0 0.09 0.1  1.0 0.0
    //     x0  y0  R    rwid vin vout
  [end]
[end]




[timestepping]
  type=be
  dt=1.0e-5
  dtmax=5.0e-2
  dtmin=5.0e-7
  endtime=1.0e2
  opts=4
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
  interval=2
[end]
