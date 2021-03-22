---
title: qpoint block
date: 2021-01-16 21:48:08
categories:
- Document
tags:
- blocks
- input file
- qpoint
---

# [qpoint] block
The block `[qpoint]` is used for the gauss point integration. This block's layout looks as follows:
```
[qpoint]
  type=gauss
  order=3
[end]
```
## options
The `type =` option specifies the type name of the gauss point integration.

`order=` specifies gauss integration order(it's not the number of gauss points, instead, it's the **order**!)


### supported gauss integration method
The full list of the available gauss integration is:
```
type=gauss
type=gausslobatto
```
