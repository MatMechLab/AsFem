*** This is an input file for the compressive neohookean model

[mesh]
  type=gmsh
  file=sample.msh
[end]



[dofs]
name=ux uy
[end]

[projection]
scalarmate=vonMises effective_plastic_strain
rank2mate=stress strain
[end]

[elmts]
  [mechanics]
    type=mechanics
    dofs=ux uy
    mate=myplastic
    domain=alldomain
  [end]
[end]

[mates]
  [myplastic]
    type=j2plasticity
    params=210.0  0.3  0.5            1.2
    //     E      nu   yield stress   hardening modulus
  [end]
[end]



[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux
    boundary=left
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dofs=uy
    boundary=left
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dofs=ux
    value=1.0*t
    boundary=right
  [end]
[end]

[output]
type=vtu
interval=5
[end]

[timestepping]
  type=be
  dt=2.0e-4
  endtime=3.0e-4
  adaptive=true
  optiters=9
  dtmax=1.0e-1
  dtmin=1.0e-12
[end]

[nonlinearsolver]
  type=nr
  maxiters=80
  r_rel_tol=5.0e-8
  r_abs_tol=4.5e-7
  solver=mumps
[end]

[job]
  type=transient
  debug=dep
[end]
