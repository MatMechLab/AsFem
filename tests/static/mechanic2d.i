// this is a test input file for mesh generation test

[mesh]
  type=asfem
  dim=2
  xmax=5.0
  ymax=5.0
  nx=100
  ny=100
  meshtype=quad9
[end]

[dofs]
name=disp_x disp_y
[end]

[qpoint]
  // for quad9 mesh, the order must>=4 !!!
  type=gauss
  order=4
[end]

[elmts]
  [elmt1]
    type=mechanics
    dofs=disp_x disp_y
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
[end]

[projection]
name=reacforce_x reacforce_y
[end]

[bcs]
  [fixbottomx]
    type=dirichlet
    dof=disp_x
    value=0.0
    boundary=bottom
  [end]
  [fixbottomy]
    type=dirichlet
    dof=disp_y
    value=0.0
    boundary=bottom
  [end]
  [loadX]
    type=dirichlet
    dof=disp_x
    value=0.2
    boundary=top
  [end]
  [loadY]
    type=dirichlet
    dof=disp_y
    value=0.1
    boundary=top
  [end]
[end]




[job]
  type=static
  debug=dep
[end]