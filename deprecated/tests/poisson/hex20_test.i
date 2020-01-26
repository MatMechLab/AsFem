*** This is a test input file for debug

[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  zmin=0.0
  zmax=5.0
  nx=20
  ny=20
  nz=20
  meshtype=hex20
[end]

[dofs]
name=phi
[end]


[elmts]
  [phi]
    type=poisson
    dofs=phi
    mate=const
    block=alldomain
  [end]
[end]

[mates]
  [const]
    type=constpoisson
    params=1.0 1.0
  [end]
[end]

[linearsolver]
  type=cg
[end]



[bcs]
  [left]
    type=dirichlet
    dof=phi
    boundary=left
    value=0.0
  [end]
  [right]
    type=dirichlet
    dof=phi
    boundary=right
    value=0.0
  [end]
  [bottom]
    type=dirichlet
    dof=phi
    boundary=bottom
    value=0.0
  [end]
  [top]
    type=dirichlet
    dof=phi
    boundary=top
    value=0.0
  [end]
  [back]
    type=dirichlet
    dof=phi
    boundary=back
    value=0.0
  [end]
  [front]
    type=dirichlet
    dof=phi
    boundary=front
    value=0.0
  [end]
[end]

[job]
  type=static
  debug=dep
[end]