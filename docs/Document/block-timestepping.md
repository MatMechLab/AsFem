---
title: timestepping block
date: 2021-01-16 21:05:08
categories:
- Document
tags:
- blocks
- input file
- timestepping
---

# [timestepping] block
The block `[timestepping]` is used for the transient analysis. This block's layout looks as follows:
```
[timestepping]
  type=be
  dt=1.0e-5
  time=1.0e3
  adaptive=true
  optiters=3
  growthfactor=1.1
  cutfactor=0.85
  dtmin=1.0e-12
  dtmax=1.0e1
[end]
```
## options
The `type =` option specifies the name of the time integration method one wants to use.

`dt=` specifies initial $\Delta_{t}$. `time=` defines the final time.

`adaptive=true` enables the adaptive time stepping. It's disabled by default.

`optiters=` determines the critical nonlinear iterations for increasing or decreasing the $\Delta_{t}$.

`growthfactor=` and `cutfactor=` represent the growth factor and cut back factor for adaptive time stepping, respectively.

`dtmin` and `dtmax` determine the minimum and maximum value of $\Delta_{t}$.

### supported time stepping method
The full list of the available time stepping method is:
```
type=be
type=cn
type=bdf2
```
