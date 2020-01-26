*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=5.0
  ymin=0.0
  ymax=5.0
  nx=80
  ny=80
  meshtype=quad8
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
    params=2.0 1.0
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
[end]

[job]
  type=static
  debug=true
[end]