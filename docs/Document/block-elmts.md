---
title: elmts block
date: 2021-01-16 10:35:08
categories:
- Document
tags:
- blocks
- input file
- elmts
---

# [elmts] block
The block `[elmts]` is used to describe the model that we plan to use according to your particular problem. This block's layout looks as follows:
```
[elmts]
  [subelmts1]
    type=elment-type2
    dofs=dof1 dof2
    mate=mate-name1
    domain=geometry-domain-name1
  [end]
  [subelmts2]
    type=elment-type2
    dofs=dof1 dof2
    mate=mate-name2
    domain=geometry-domain-name2
  [end]
  ...
[end]
```
## options
The `type =` option specifies the name of the element type (or the physical model) one wants to use.

`dofs=` specifies which DoFs we want to use, it should be noted that one of the names in your `[dofs]` block must be the name we used here.

`mate=` gives the material name, which should be the **block name** in your material block, **not** the material type name!!! This option can be ignored, then your element will not call any material calculation.

`domain=` determines which domain will be applied to the current element. Users don't usually need this option, then all of your mesh domains will be used by default.

### supported element type
The full list of the available element type is:
```
type=poisson
type=mechanics
type=diffusion
type=cahnhilliard
type=user1[,...,user20]
```
