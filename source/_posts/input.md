---
title: Input file
mathjax: false
date: 2020-04-16 16:30:08
categories:
- Documents
tags:
- Documents
- Input
---

# Before we start
Currently, some available demos can be found here:

[Demo](https://space.bilibili.com/100272198/channel/detail?cid=90241)

## Write a simple input file:

structure of the block in your input file should look like:
```
[block_name]
type=...
option1=...
option2=...
[end]
```
Here we use the [blockname]/[end] bracket pair to set the properties we want. Thus each block pair must end up with an '[end]'. Remember, even there are sub-blocks inside one big block, each block itself still must end up with an '[end]'. Otherwise, your input file will complain about errors to you!

For one minimal input file, you need:
```
[mesh]
...
[end]
```

and the name or name list of the defined degrees of freedom (DoFs), which will be displayed in the Paraview when you open the result file(vtu file). Besides, the name of related DoFs is required when you want to apply the boundary condition.
```
[dofs]
...
[end]
```

and the element or the module you want to use:
```
[elmts]
...
[end]
```

and the information for the analysis(we call it the [job] block)
```
[job]
...
[end]
```

In summary, you need at least:
```
[mesh]
...
[end]
[dofs]
...
[end]
[elmts]
...
[end]
[job]
...
[end]
```
