[mesh]
  type=asfem
  dim=2
  nx=5
  ny=5
  meshtype=quad4
[end]

[dofs]
name=c1 c2 c3
[end]

[elmts]
  [myelmt]
    type=user2
    dofs=c1 c2 c3
  [end]
[end]


[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-9
  r_abs_tol=4.6e-7
  //solver=superlu
[end]



[bcs]
  [fixc1]
    type=dirichlet
    dof=c1
    value=0.0
    boundary=left
  [end]
  [fixc1r]
    type=dirichlet
    dof=c1
    value=1.0
    boundary=right
  [end]
  [fixc2]
    type=dirichlet
    dof=c2
    value=0.0
    boundary=left
  [end]
  [fixc2r]
    type=dirichlet
    dof=c2
    value=1.0
    boundary=right
  [end]
  [fixc3]
    type=dirichlet
    dof=c3
    value=0.0
    boundary=left
  [end]
  [fixc3r]
    type=dirichlet
    dof=c3
    value=1.0
    boundary=right
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
