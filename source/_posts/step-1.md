---
title: Step 1 - Mesh generation
mathjax: true
date: 2021-01-02 17:30:08
categories:
- Tutorial
tags:
- tutorial
- input file
- mesh
---

# Introduction

Before we begin the actual FEM simulation, we must define our computation domain and discretize it into multiple subdomains. For this purpose, the `[mesh]` block is introduced.

# 1D example

Let's take the 1D solid line as an example, this line can be discretized into several 1D Lagrange mesh as follows:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=10
  meshtype=edge2
[end]
```
## Options
The `type= `  option specifies the type of mesh generation we plan to use. We are offering two kinds of mesh generation in AsFem. The first one is the built-in mesh generation for a regular domain, i.e. the straight line (1d), the rectangle domain (2d), and the cubic domain (3d). For the second one, users can import their favorite mesh from other packages, like gmsh ( `type=gmsh` ) or netgen ( `type=gmsh2` ).

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

## Complete mesh block
Then the complete `[mesh]` block should look like:
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

# 2D and 3D examples

Similarly, for 2D and 3D cases, one can use:
```
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=10
  ny=10
  meshtype=quad4
[end]
```
and
```
[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=1.0
  nx=10
  ny=10
  nz=10
  meshtype=hex8
[end]
```
Or, one can also use:
```
[mesh]
  type=asfem
  dim=3
  nx=10
  ny=10
  nz=10
  meshtype=hex8
[end]
```
then the unit `[0,1]x[0,1]x[0,1]` 3D domain will be used by default.

# First try in your AsFem

Now, lets try your first example in AsFem. You can create a new txt file or simply run the following commands:
```
nano firstrun.i
```
or
```
vim firstrun.i
```
then copy and paste the following `[mesh]` block into your `firstrun.i`:
```
[mesh]
  type=asfem
  dim=3
  nx=10
  ny=10
  nz=10
  meshtype=hex8
  savemesh=true
[end]
```
save it and then execute your `AsFem` as follows:
```
asfem -i firstrun.i
```
or in parallel:
```
mpirun -np 4 asfem -i firstrun.i
```


The complete input files can be fund in `examples/tutorial`.
