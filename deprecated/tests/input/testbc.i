*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=2.0
  ymin=0.0
  ymax=2.0
  nx=20
  ny=20
  meshtype=quad4
[end]

[dofs]
name=u v w
[end]

[elmts]
	[elmt1sss]
	  type=solid2d
	  dofs=u v w
	[end]

	[elmt2]
	  type=solid2d
	  dofs=u v w
	  domain=test
	  //mate=user0
	[end]

	[elmt3]
	  type=user1
	  dofs=u v
	  //mate=user1
	[end]
[end]

[bcs]
  [bc1]
    type=dirichlet
    dof=u
    boundary=left
  [end]

  [bc2]
    type=user15
    dof=w
    boundary=top
  [end]
  [bc3]
    type=neumann
    dof=v
    boundary=right
  [end]
[end]