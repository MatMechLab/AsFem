[mesh]
  type=gmsh
  file=sample.msh
  savemesh=true
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

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-9
  r_abs_tol=4.6e-7
  solver=superlu
[end]



[bcs]
  [fixux]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=bottom
  [end]
  [fixuy]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom
  [end]
  [load]
    type=dirichlet
    dof=uy
    value=0.1
    boundary=top
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
