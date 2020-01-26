*** 2d linear elastic problem

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  nx=20
  ny=20
  meshtype=quad4
[end]

[nonlinearsolver]
  //type=nrregulalinesearch
  //type=nrbisectlinesearch
  type=nr
  stol=1.0e-4
  r_abs_tol=2.5e-7
  e_abs_tol=3.5e-15
  maxiters=500
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [solids]
    type=mechanics
	  dofs=ux uy
    mate=sv
    block=alldomain
  [end]
[end]

[mates]
  [sv]
    type=saintvenant
    params=1.0e5 0.25
    // params= E nu(Youngs modulus and poisson ratio)
  [end]
[end]


[bcs]
  [fix_ux]
    type=preset
    dof=ux
    value=0.0
    boundary=left
  [end]
  [fix_uy]
    type=preset
    dof=uy
    value=0.0
    boundary=bottom
  [end]
  [loadux]
    type=dirichlet
    dof=ux
    value=0.1
    boundary=right
  [end]
  [shearux]
    type=dirichlet
    dof=uy
    value=0.15
    boundary=top
  [end]
[end]

[job]
  type=static
  debug=dep
  projection=true
[end]