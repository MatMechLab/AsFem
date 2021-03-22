

# Introduction
In [step-1](https://yangbai90.github.io/AsFem/2021/01/02/step-1/), our mesh was defined. Some necessary information, however, is still missing for a simple FEM analysis. Therefore, we will try to define our degrees of freedom (DoFs) in this step and also the element of our model. In the end, to obtain the final solution, we will apply the related boundary conditions.

# The poisson equation
The issue we want to solve is the Poisson linear equation that reads as follows:
$$
\begin{equation}
k\nabla^{2}\phi=F
\label{eq:poisson}
\tag{1}
\end{equation}
$$
where $k$ and $F$ denote the model's coefficients. Below are the related boundary conditions:
$$
\begin{equation}
k\nabla\phi\cdot\vec{n}=0\qquad\mathrm{on}\quad\partial\Omega_{N}
\label{eq:neumann}
\tag{2}
\end{equation}
$$
and
$$
\begin{equation}
\phi=\phi_{g}\qquad\mathrm{on}\quad\partial\Omega_{D}
\label{eq:dirichlet}
\tag{3}
\end{equation}
$$
where the *Dirichlet* boundary condition and *Neumann* boundary condition are defined by the subscripts $D$ and $N$.


# Define the degree of freedom (DoF)

The degree of freedom (DoF) or the degrees of freedom (DoFs) can be used to define the name of each DoF and also to apply the necessary boundary conditions (`[bcs]`), elements (`[elmts]`), and so on. The `[dofs]` block looks like below:
```
[dofs]
name=dof1 dof2 dof3 ...
[end]
```
## Options
The `name=`  option specifies the name of each DoF. One should keep in mind that, the order of the name indicates the index of each DoFs. For instance, we need two displacements, namely `disp_x` and `disp_y`, if we want to do a 2D elastic analysis. The block of `[dofs]` should therefore be specified as:
```
[dofs]
name=disp_x disp_y
[end]
```
where `disp_x` is the first DoF(index=1), `disp_y` is the second DoF(index=2). That's all, `name=` is the only option in `[dofs]` block, nothing else.

For the Poisson equation, because there is only one DoF involved, the final expression of `[dofs]` should be:
```
[dofs]
name=phi
[end]
```

# Element for Poisson equation
The DoF is ready now, but the model in Eq.$\eqref{eq:poisson}$ is still missing. Thereby, we introduce the `[elmts]` block for this purpose. This block looks like below:
```
[elmts]
  [elmt1]
    type=poisson
    dofs=phi
    mate=mymate
  [end]
[end]
```
where `type=` option specifies the element we want to use, it could be either the built-in elements of AsFem or the user-defined-element (`UEL`). The DoFs that will be used in this element are defined by `dofs=`. `mate=` gives the name of the material block that we want to use. Once the `[elmts]` block is given, the model we defined in Eq.$\eqref{eq:poisson}$ is ready.

## Material properties
For the coefficients $k$ and $F$, namely the material properties, they can be calculated or defined via the `[mates]` block as follows:
```
[mates]
  [mymate]
    type=constpoisson
    params=1.0 1.0e1
  [end]
[end]
```
where `type=` specifies the material type name defined in AsFem. `params=` defines the parameters we want to use in our model, in this case, $k=1.0$ and $F=10.0$ will be used.

# Boundary conditions
The boundary conditions, as mentioned in Eq. $\eqref{eq:dirichlet}$ and $\eqref{eq:neumann}$, can be applied via the `[bcs]` block. In our case, the *Neumann* boundary condition in Eq.$\eqref{eq:neumann}$ is zero, therefore, only the *Dirichlet* boundary condition need to be considered:
```
[bcs]
  [fixleft]
    type=dirichlet
    dof=phi
    value=0.1
    boundary=left
  [end]
  [fixright]
    type=dirichlet
    dof=phi
    value=0.5
    boundary=right
  [end]
[end]
```
where `type=` specifies the different types of boundary conditions supported by AsFem. `dof=` denotes the *name* of DoF we want to apply the given boundary conditions. In our case, we constrain the value of $\phi$ on the *left* and *right* side of a rectangle domain to be 0.1 and 0.5, respectively.


# Static analysis
Until now, all the model and boundary conditions are ready. To start the FEM calculation, we need a `[job]` block to tell AsFem which kind of analysis we want. For the static analysis in this case, it can be given as follows:
```
[job]
  type=static
[end]
```
if one wants to see how the iteration information changes, one can use:
```
[job]
  type=static
  debug=dep
[end]
```
where `debug=` option enables some basic information output in your terminal. If you don't want to see too many outputs, then you can use `debug=false`.



# Run it in AsFem
Now, let's try your second example in AsFem. You can create a new text file and name it as step2.i or whatever you like. Then copy the following lines into your input file:
```
[mesh]
  type=asfem
  dim=2
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=phi
[end]

[elmts]
  [elmt1]
    type=poisson
    dofs=phi
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=constpoisson
    params=1.0 1.0e1
  [end]
[end]

[bcs]
  [fixleft]
    type=dirichlet
    dof=phi
    value=0.1
    boundary=left
  [end]
  [fixright]
    type=dirichlet
    dof=phi
    value=0.5
    boundary=right
  [end]
[end]


[job]
  type=static
[end]
```
You can also find the complete input file in `examples/tutorial`.


If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):
![](step2.jpeg)



Wait a minute, what should I do if I want to solve a Poisson equation in 3D?


The answer is ...... quite simple, just change your mesh to 3D like this(the complete input file is `step2-3d.i`):
```
[mesh]
  type=asfem
  dim=3
  nx=50
  ny=50
  nz=50
[end]
```
then you will see:
![](step2-3d.jpeg)
