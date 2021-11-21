[mesh]
  type=gmsh
  file=rect.msh
  //file=sample.msh
[end]


[qpoint]
  type=gauss
  order=3
[end]

[dofs]
name=u
[end]

[elmts]
  [poisson]
    type=poisson
    dofs=u
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constpoisson
    params=1.0 1.0
    //     
  [end]
[end]

[nonlinearsolver]
  type=nr
  maxiters=50
  r_rel_tol=5.0e-10
  r_abs_tol=2.5e-7
  //solver=superlu
[end]



[output]
  type=vtu
  interval=1
[end]


[bcs]
  [fixux]
    type=dirichlet
    dofs=u
    value=0.0
    boundary=left right bottom top
  [end]
[end]



[postprocess]
  [top]
    type=area
    side=top
  [end]
  [left]
    type=area
    side=left
  [end]
  [area]
    type=volume
    domain=block
  [end]
[end]

[job]
  type=static
  debug=dep
[end]
