// this is a test input file for mesh generation test

[mesh]
  type=abaqus
  //file=Quad4-MultipleSurf.inp
  //file=quad8.inp
  //file=tet4.inp
  //file=tet10.inp
  //file=tri3.inp
  file=Tri6-SurfaceSet.inp
  savemesh=true
[end]

[dofs]
name=c
[end]

[elmts]
  [elmt1]
    type=diffusion
    dofs=c
    mate=mate1
  [end]
[end]

[mates]
  [mate1]
    type=constdiffusion
    params=1.0e2
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  time=2.0e-5
  adaptive=true
  optiters=3
[end]

[projection]
vectormate=gradc
[end]

//[ics]
//  [randc]
//    type=random
//    dof=c
//    params=1.0 2.0
//  [end]
//[end]

[bcs]
  //[fix]
  //  type=nodaldirichlet
  //  dofs=c
  //  value=1.0
  //  boundary=Node-Left
  //[end]
  [flux]
    type=nodalneumann
    dofs=c
    value=-0.1
    boundary=Node-Left
  [end]
[end]

[job]
  type=transient
  debug=dep
[end]
