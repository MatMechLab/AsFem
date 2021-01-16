---
title: timestepping block
date: 2021-01-16 21:15:08
categories:
- Document
tags:
- blocks
- input file
- ics
---

# [ics] block
The block `[ics]` is used to apply the initial condition for the DoFs. This block's layout looks as follows:
```
[ics]
  [ic1]
    type=const
    dof=dof1
    domain=domain-name
    params=val1 val2 ...
  [end]
  ...
[end]
```
## options
The `type =` option specifies the name of the initial condition type one wants to use.

`dof=` specifies the name of the DoF.

`domain=` indicates the domain name where the initial condition will be applied.

`params=` gives the parameters which will be used in the initial condition calculation.



### supported initial condition type
The full list of the available initial condition typ is:
```
type=const
type=random
type=circle
type=sphere
type=rectangle
type=cubic
```
