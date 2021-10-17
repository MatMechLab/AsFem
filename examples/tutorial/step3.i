[mesh]
  type=asfem
  dim=2
  xmax=5.0
  ymax=5.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
  [end]
[end]

[bcs]
  [fixbottomX]
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=bottom
  [end]
  [fixbottomY]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=bottom
  [end]
  [loadY]
    type=dirichlet
    dofs=uy
    value=0.1
    boundary=top
  [end]
[end]

[projection]
name=reacforce_x reacforce_y
scalarmate=vonMises
rank2mate=stress
[end]

[job]
  type=static
  debug=dep
[end]