*** This is an input file for tesing

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=2.0
  ymin=0.0
  ymax=2.0
  nx=2
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
	  mate=mate1
	[end]

	[elmt3]
	  type=user1
	  dofs=u v
	  mate=user1
	[end]
[end]

[mates]
  [mate1]
    type=const
  [end]
  [user1]
    type=freeenergy
  [end]
[end]

[ics]
  [ic1]
    type=const
    dof=u
    block=all

  [end]

  [ic2]
    type=random
    dof=v
    block=1
  [end]

  [ic3]
    type=user14
    dof=w
    block=all1
  [end]
[end]
