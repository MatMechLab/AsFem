---
title: nonlinear block
date: 2021-10-29 21:28:08
categories:
- Document
tags:
- blocks
- input file
- projection
---

# [projection] block
The block `[projection]` is used for projecting the quantities of the gauss point to the nodal point. This block's layout looks as follows:
```
[projection]
  name=x1 x2 x3
  scalarmate=mate1 mate2 ...
  vectormate=vmat1 vmat2 ...
  rank2mate=stress strain ...
  rank4mate=jacobian ...
[end]
```
## options
The `name=` option specifies the name of the scalar value you want to project, it should be calculated/defined in your element, it is independent with the material system.

`scalarmate=` specifies the name of scalar type materials you want to project, which should be defined/calculated in your materials. For example, `Mate.ScalarMaterials("myscalar")=0.0`, where "myscalar" is the name for the vector type material.

`vectormate=` specifies the name of vector type materials you want to project, which should be defined/calculated in your materials. For example, `Mate.VectorMaterials("myvector")=1.0`, where "myvector" is the name for the vector type material.

`rank2mate=` specifies the name of rank-2 tensor type materials you want to project, which should be defined/calculated in your materials. For example, `Mate.Rank2Materials("mystress")=1.0`, where "mystress" is the name for the rank-2 tensor type material.


`rank4mate=` specifies the name of rank-4 tensor type materials you want to project, which should be defined/calculated in your materials. For example, `Mate.Rank4Materials("myjacobian")=1.0`, where "myjacobian" is the name for the rank-4 tensor type material.


Once you put the correct name of the material properties there, AsFem will automatically save them into the result file. Then you can check them easily in the Paraview.
