## Write a simple input file:

In your input file, the layout of the block should look like:
```
[block_name]
type=...
option1=...
option2=...
[end]
```
here, to set the properties we want, we use the `[blockname]/[end]` bracket pair. Every block pair must therefore end up with an `[end]`. Note, even inside one wide block there are sub-blocks, each (sub)block itself always has to end up with an `[end]`. Otherwise, the input file will complain to you about mistakes!


For one minimal input file, you need:
```
[mesh]
...
[end]
```

and also the *Specified Degrees of Freedom* (**DoFs**) name or name list that will be displayed in the Paraview when you open the result file (vtu file). Furthermore, when you want to apply the boundary condition, the name of related **DoFs** is required.
```
[dofs]
...
[end]
```

and the element or module you'd like to use:
```
[elmts]
...
[end]
```

as well as the necessary information for the analysis(we call it the `[job]` block)
```
[job]
...
[end]
```

To sum up, you at least need:
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
However, one can also use the `--read-only` option to force AsFem to only read the input file but ignore the FEM analysis workflow.

The list of all the blocks in full is:
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

[mates]
...
[end]

[bcs]
...
[end]

[ics]
...
[end]

[projection]
...
[end]

[timestepping]
...
[end]

[nonlinearsolver]
...
[end]

[qpoint]
...
[end]

[job]
...
[end]
```

In the following steps, we will introduce you the details of each block and their functions.

All the input file can be found in the `examples/tutorial` folder.
