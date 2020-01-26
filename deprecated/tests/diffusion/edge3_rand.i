*** simple diffusion for 1d case

[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=2.0
  nx=500
  meshtype=edge3
[end]

[dofs]
name=conc
[end]

[elmts]
  [phi]
    type=diffusion
    dofs=conc
    mate=const
    block=alldomain
  [end]
[end]

[mates]
  [const]
    type=constdiffusion
    params=2.5
  [end]
[end]

[ics]
  [rand]
    type=random
    dof=conc
    values=0.10 0.90
    block=alldomain
  [end]
[end]




[job]
  type=transient
  debug=true
  dt=1.0e-5
  endtime=5.0e-1
  interval=2
[end]