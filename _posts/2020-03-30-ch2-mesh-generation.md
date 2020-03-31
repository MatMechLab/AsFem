---
layout: post
title: Chapter 2-Mesh generation
permalink: /documents/
---
For the mesh generation, the [mesh] block is required.
For example, if one wants to generate the 1d Lagrange mesh, you can use:
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
if one wants to save the generated mesh, then the option 'savemesh=true' is required:
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
The mesh will be saved as a .vtu file, which should have the name: your_input_file_name+'_mesh.vtu'(.i is removed from your input file name). For example, if your input file is: test.i, then the mesh file name is: test_mesh.vtu. 

Additionally, AsFem will also print this message in your console:
```
**************************************************************
*** Start to read input file ...                           ***
***   start to crate mesh ...                              ***
***     save mesh to                   hex20m1_mesh.vtu    ***
***   mesh generation finished !                           ***
```

Moreover, if one wants to print out the mesh information in the terminal, one can use the 'printmesh=true' option. For example, if we want to print the information of a 1D edge3 type Lagrange mesh, we can use:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=5
  meshtype=edge4
  savemesh=false
  printmesh=true
[end]
```
then the output should look like:
```
***-------------------------------------------------------------------***
*** Mesh information summary:                                         ***
***  Nodes=       16, Elmts=        7, NodesPerBulkElmt=  4           ***
***  Max dim= 1, Min dim= 0, PhyGroup=    3, Meshtype= edge4, Order=3 ***
***  Physical id                    Phsical Name             Elmts    ***
***        1                                left                 1    ***
***        2                               right                 1    ***
***        3                           alldomain                 5    ***
***-------------------------------------------------------------------***
```
if one wants to print out the details of his mesh, then 'printmesh=dep' option should be used. Consequently, the output should look like:
```
***-------------------------------------------------------------------***
*** Mesh information summary:                                         ***
***  Nodes=       16, Elmts=        7, NodesPerBulkElmt=  4           ***
***  Max dim= 1, Min dim= 0, PhyGroup=    3, Meshtype= edge4, Order=3 ***
***  Physical id                    Phsical Name             Elmts    ***
***        1                                left                 1    ***
***        2                               right                 1    ***
***        3                           alldomain                 5    ***
***  Physical group information (ID and element ID)                   ***
***  phyname=                     left, element id=        1          ***
***  phyname=                    right, element id=        2          ***
***  phyname=                alldomain, element id=        3          ***
***  phyname=                alldomain, element id=        4          ***
***  phyname=                alldomain, element id=        5          ***
***  phyname=                alldomain, element id=        6          ***
***  phyname=                alldomain, element id=        7          ***
***  Element connectivity information(element id: node index):        ***
***  elmt id=        1:       1                                       ***
***  elmt id=        2:      16                                       ***
***  elmt id=        3:       1        2        3        4            ***
***  elmt id=        4:       4        5        6        7            ***
***  elmt id=        5:       7        8        9       10            ***
***  elmt id=        6:      10       11       12       13            ***
***  elmt id=        7:      13       14       15       16            ***
***  Node coornidates (node id, x, y, z and weight)                   ***
***          1:   0.0000e+00,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          2:   6.6667e-02,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          3:   1.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          4:   2.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          5:   2.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          6:   3.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          7:   4.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          8:   4.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          9:   5.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         10:   6.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         11:   6.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         12:   7.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         13:   8.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         14:   8.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         15:   9.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         16:   1.0000e+00,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***-------------------------------------------------------------------***
```
Since **printmesh=dep** option will print out all the information, includes the nodal coordinates and elemental connectivity, it is suggested **not** use this option in your [mesh] block. If one wants to check the mesh, 'savemesh=true' is already good enough!

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
then the unit([0,1] in 1D, [0,1]x[0,1] in 2D, [0,1]x[0,1]x[0,1] in 3D) domain will be used by default.

Currently, AsFem supports the following kinds of built-in mesh:

1D--> edge2, edge3, edge4

2D--> quad4, quad8, quad9

3D--> hex8, hex20, hex27