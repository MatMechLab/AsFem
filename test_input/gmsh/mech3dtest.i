[mesh]
  type=asfem
  dim=3
  xmax=1
  ymax=2
  zmax=3
  nx=10
  ny=20
  nz=30
  meshtype=hex8
[end]

[dofs]
name=ux uy uz
[end]

[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy uz
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
  [end]
[end]

[nonlinearsolver]
  //type=helper
  type=nr
  maxiters=50
  r_rel_tol=1.0e-12
  r_abs_tol=2.6e-7
[end]



[bcs]
  [fixux]
    //type=helper
    type=dirichlet
    dofs=ux
    value=0.0
    boundary=back
  [end]
  [fixuy]
    type=dirichlet
    dofs=uy
    value=0.0
    boundary=back
  [end]
  [fixuz]
    type=dirichlet
    dofs=uz
    value=0.0
    boundary=back
  [end]
  [load]
    type=dirichlet
    dofs=uz
    value=0.1
    boundary=front
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]


[postprocess]
  [arealeft]
    //type=helper
    type=area
    side=left
  [end]
  [areabottom]
    type=area
    side=bottom
  [end]
  [areaback]
    type=area
    side=back
  [end]
  [volume]
    type=volume
    domain=alldomain
  [end]
  [nodeux]
    type=nodevalue
    dof=ux
    nodeid=1211
  [end]
  [nodeuy]
    type=nodevalue
    dof=uy
    nodeid=1211
  [end]
  [nodeuz]
    type=nodevalue
    dof=uz
    nodeid=1211
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
