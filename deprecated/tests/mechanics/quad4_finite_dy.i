*** 2d linear elastic problem

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [solids]
    type=mechanics
    dofs=ux uy
    mate=finite
    block=alldomain
  [end]
[end]

[mates]
  [finite]
    type=finitestrain
    params=1.0e2 0.25 0
    // params= E nu(Youngs modulus and poisson ratio)
  [end]
[end]


[bcs]
  [fix_ux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=left
  [end]
  [fix_uy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom
  [end]
  [loadux]
    type=dirichlet
    dof=ux
    value=0.2*t
    boundary=right
  [end]
  [shearux]
    type=dirichlet
    dof=uy
    value=0.5*t
    boundary=top
  [end]
[end]

[job]
  type=transient
  debug=true
  dt=1.0e-3
  dtmax=1.0e1
  adaptive=true
  projection=true
[end]