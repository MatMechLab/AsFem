*** This is an input file for the compressive neohookean model

[mesh]
  type=gmsh
  file=logo.msh
[end]



[dofs]
name=ux uy
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]

[elmts]
  [mechanics1]
    type=mechanics
    dofs=ux uy
    mate=smallstrain1
    domain=matrix1
  [end]
  [mechanics2]
    type=mechanics
    dofs=ux uy
    mate=smallstrain2
    domain=matrix2
  [end]
  [mechanics3]
    type=mechanics
    dofs=ux uy
    mate=smallstrain3
    domain=inclusion
  [end]
[end]

[mates]
  [smallstrain1]
    type=incresmallstrain
    params=10.0 0.3
  [end]
  [smallstrain2]
    type=incresmallstrain
    params=20.0 0.3
  [end]
  [smallstrain3]
    type=incresmallstrain
    params=500.0 0.3
  [end]
[end]



[bcs]
  [FixUx]
    type=dirichlet
    dof=ux
    boundary=left
    value=0.0
  [end]
  [FixUy]
    type=dirichlet
    dof=uy
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=ux
    value=1.0*t
    boundary=right
  [end]
  [loadUy]
    type=dirichlet
    dof=uy
    value=2.0*t
    boundary=top
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-3
  endtime=1.0e-2
  adaptive=true
  optiters=3
  dtmax=1.0e-1
[end]

[job]
  type=transient
  debug=dep
[end]
