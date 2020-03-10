*** This is an input file for 2d thermal analysis (Joule heat from current)
[mesh]
  type=gmsh
  file=fiber.msh
[end]

[dofs]
name=T
[end]



[elmts]
  [thermal1]
    type=thermal
    dofs=T
    mate=mate1
    domain=inner
  [end]
  [thermal2]
    type=thermal
    dofs=T
    mate=mate2
    domain=outer
  [end]
[end]

[mates]
  [mate1]
    type=currentthermal
    params=1.0e-3   5.0                 5.0e-2      1.0e-1  1500.0
    //     Capacity thermal-diffusivity resistivity current density
  [end]
  [mate2]
    type=currentthermal
    params=4.0e-3   12.5                10.0e-2     1.0e-1  2000.0
    //     Capacity thermal-diffusivity resistivity current density
  [end]
[end]

[ics]
  [T0]
    type=const
    dof=T
    params=25.0
  [end]
[end]




[timestepping]
  type=be
  dt=1.0e-5
  dtmax=2.0e1
  dtmin=5.0e-7
  endtime=1.0e4
  adaptive=true
[end]
[nonlinearsolver]
  type=newtonls
  r_abs_tol=5.0e-7
  r_rel_tol=1.5e-8
  maxiters=15
[end]


[job]
  type=transient
  debug=dep
[end]
