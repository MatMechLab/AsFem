[mesh]
  type=asfem
  dim=2
  xmax=10.0
  ymax=10.0
  nx=100
  ny=100
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
    type=user1
    params=210.0 0.3 0.2
    //     E0    nu  delta
  [end]
[end]

[bcs]
  [fix]
    type=dirichlet
    dofs=ux uy
    value=0.0
    boundary=bottom
  [end]
  [loading]
    type=dirichlet
    dofs=uy
    value=0.02
    boundary=top
  [end]
[end]

[projection]
scalarmate=MyE vonMises
[end]


[job]
  type=static
  debug=dep
[end]