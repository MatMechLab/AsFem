*** user-defined element for nonlinear geometric problem (elastic)

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  nx=100
  ny=100
  meshtype=quad9
[end]

[dofs]
name=ux uy
[end]

[nonlinearsolver]
  //type=nrbisectlinesearch
  //type=nrregulalinesearch
  type=nr
  //stol=1.0e-4
[end]

[linearsolver]
  type=lu
[end]

[elmts]
  [solids]
    type=user1
	  dofs=ux uy
    mate=umat1
    block=alldomain
  [end]
[end]

[mates]
  [umat1]
    type=user1
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