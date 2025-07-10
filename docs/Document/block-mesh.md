---
title: mesh block
date: 2021-01-15 16:03:08
categories:
- Document
tags:
- blocks
- input file
- mesh
---

# [mesh] block
The format of the block is:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=10
  meshtype=edge2
  savemesh=true
[end]
```
## [mesh] block options
`type= `  option specifies the type of mesh generation we plan to use. We are offering two kinds of mesh generation in AsFem. The first one is the built-in mesh generation for a regular domain, i.e. the straight line (1d), the rectangle domain (2d), and the cubic domain (3d). For the second one, users can import their favorite mesh from other packages, like gmsh ( `type=gmsh` ) or netgen ( `type=gmsh2` ).

`dim=`  determines the domain's dimension, which should be 1, 2, or 3.

`xmin=` and `xmax=` denote the size of the domain, you will need `ymin=` and `ymax=` for the 2D case, and `zmin=` and `zmax=` in the 3D case. One can also **ignore** these options, by default, the size of the domain will be unit, namely `[1]` in 1D, `[1,1]` in 2D, `[1,1,1]` in 3D.

Simultaneously, `nx`, `ny`, and `nz` represent the number of mesh along these three axes, respectively.

`meshtype=` option offers the choices of different kinds of mesh, for instance, the second order Lagrange mesh in 1D case can be obtained via `meshtype=edge3`. Currently, AsFem offers:
```
edge2,edge3,edge4 // in 1D case
quad4,quad8,quad9 // in 2D case
hex8, hex20,hex27 // in 3D case
```

If one want to save the created mesh, one will need the `savemesh=true` option. The mesh will be saved as a *.vtu* file, which should be named as 'your_input_file_name'+'_mesh.vtu' (*.i* is removed from your input file name). For example, if your input file is: *test.i*, then the mesh file name is: *test_mesh.vtu*.
