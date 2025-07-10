---
title: dofs block
date: 2021-01-16 10:26:08
categories:
- Document
tags:
- blocks
- input file
- dofs
---

# [dofs] block
The degree of freedom (DoF) or the degrees of freedom (DoFs) can be used to define the name of each DoF and also to apply the necessary boundary conditions (`[bcs]`), elements (`[elmts]`), and so on. The `[dofs]` block looks like below:
```
[dofs]
name=dof1 dof2 dof3 ...
[end]
```
## [dofs] block option
The `name=`  option specifies the name of each DoF. One should keep in mind that, the order of the name indicates the index of each DoFs. For instance, we need two displacements, namely `disp_x` and `disp_y`, if we want to do a 2D elastic analysis. The block of `[dofs]` should therefore be specified as:
```
[dofs]
name=disp_x disp_y
[end]
```
where `disp_x` is the first DoF(index=1), `disp_y` is the second DoF(index=2).

That's all, `name=` is the only option in `[dofs]` block, nothing else.
