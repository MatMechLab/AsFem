// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=3
  xmax=2.0
  ymax=2.0
  zmax=2.0
  nx=5
  ny=5
  nz=5
  meshtype=hex8
[end]

[dofs]
name=disp_x disp_y disp_z
[end]

[qpoint]
  type=gauss
  order=2
[end]

[elmts]
  [elmt1]
    type=mechanics
    dofs=disp_x disp_y disp_z
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=linearelastic
    params=120.0 0.3
    //     E     nu
  [end]
[end]



[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=1.0e-10
  r_abs_tol=1.0e-8
  solver=mumps
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
rank4mate=jacobian
[end]

[bcs]
  [fixbottomx]
    type=dirichlet
    dofs=disp_x
    value=0.0
    boundary=bottom
  [end]
  [fixbottomy]
    type=dirichlet
    dofs=disp_y
    value=0.0
    boundary=bottom
  [end]
  [fixbottomz]
    type=dirichlet
    dofs=disp_z
    value=0.0
    boundary=bottom
  [end]
  [loadX]
    type=dirichlet
    dofs=disp_x
    value=0.2
    boundary=top
  [end]
  [loadY]
    type=dirichlet
    dofs=disp_y
    value=0.1
    boundary=top
  [end]
[end]




[job]
  type=static
  debug=dep
[end]