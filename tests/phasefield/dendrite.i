*** This is an input file for 2d phase field fracture model

[mesh]
  type=asfem
  dim=2
  xmin= 0.0
  xmax=10.0
  ymin= 0.0
  ymax=10.0
  nx=150
  ny=150
  meshtype=quad4
[end]

[dofs]
name=F dFdphi K dKdtheta T_x T_y
[end]

[projection]
  name=vonMises Hydro sigxx sigyy sigxy
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
    params=3333.33 0.01 0.04  6 90.0   1.0     1.8
    //     L       eps  delta J theta0 Conduct eta
  [end]
[end]

[ics]
  [const]
    type=const
    dof=d
    params=0.0
  [end]
[end]




[timestepping]
  type=be
  dt=1.0e-5
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
