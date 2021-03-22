---
title: bcs block
date: 2021-01-16 16:15:08
categories:
- Document
tags:
- blocks
- input file
- bcs
---

# [bcs] block
The block `[bcs]` is used to apply the different types of boundary conditions, i.e. *Dirichlet* boundary condition, *Neumann* boundary condition, *Robin* boundary condition, as well as *User-Defined-BC* (**UBC**). This block's layout looks as follows:
```
[bcs]
  [mybc1]
    type=dirichlet
    dof=dof1
    boundary=left right ...
    value=bcvalue
  [end]
[end]
```

## options
The `type =` option specifies the name of the boundary condition type one wants to use.

`dof=` specifies which DoF we want to use. It should be noted that, only one DoF can be accepted, instead of multiple DoFs used in other blocks like `[elmts]`.

`boundary=` specifies the name of the boundary, which we want to apply the related boundary condition. The name of the boundary should be defined in your mesh file. For the built-in mesh, we use `left` and `right` for the point/line/surface at `xmin` and `xmax`. Similarly, `bottom` and `top` are used for the line/surface at `ymin` and `ymax`. Then `back` and `front` are used for the surface at `zmin` and `zmax`.

`value=` specifies the boundary value we want to use. It should be a single value, instead of several numbers. If one wants to apply the time dependent boundary condition, then he can use `value=t*2.0`. Thus the boundary value will change overtime.


### supported boundary condition type
The full list of the available boundary condition type is:
```
type=dirichlet
type=neumann
```
