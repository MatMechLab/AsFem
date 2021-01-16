---
title: mates block
date: 2021-01-16 16:15:08
categories:
- Document
tags:
- blocks
- input file
- mates
---

# [mates] block
The block `[mates]` is used to calculate the material properties we plan to use according to our model, therefore, the material properties will be used by the `[elmts]` block. This block's layout looks as follows:
```
[mates]
  [mymate1]
    type=linearelastic
    params=210.0 0.3
  [end]
  [mymate2]
    type=linearelastic
    params=410.0 0.2
  [end]
  ...
[end]
```
It should be noted that, the name of the block `mymate1` should be the same as the one you give in `[elmts]` subblock(the `mate=mymate1` option).

## options
The `type =` option specifies the name of the element type (or the physical model) one wants to use.

`params=` specifies the parameters we want to use in the material property calculation.


### supported material type
The full list of the available element type is:
```
type=constpoisson
type=linearelastic
type=cahnhilliard
```
