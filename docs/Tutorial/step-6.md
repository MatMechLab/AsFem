---
title: Step 6 - Solve the linear momentum balance equation
mathjax: true
date: 2021-10-29 20:20:08
categories:
- Tutorial
tags:
- tutorial
- input file
- mesh
- dofs
- elmts
- mates
- stress
- linear momentum balance
- timestepping
---

# Introduction
In this step, we will try to solve the [linear momentum balance](https://en.wikiversity.org/wiki/Continuum_mechanics/Balance_of_linear_momentum) equation for the linear elastic problem.

# The linear momentum balance equation
The equation can be read as follows:
$$
\begin{equation}
\rho\frac{\partial v}{\partial t}=\mathbf{\nabla}\cdot\mathbf{\sigma}+\rho\mathbf{b}
\label{eq:stress-equation}
\tag{1}
\end{equation}
$$
with $\rho$ being the density, $\mathbf{b}$ denoting the body force (i.e., gravity), $\mathbf{\sigma}$ being the Cauchy stress tensor. The boundary condition for this problem is defined as follows:
$$
\begin{equation}
\mathbf{u}=u_{g}\qquad\mathrm{on}\quad\partial\Omega_{D}\qquad\qquad
\mathbf{\sigma}\cdot\vec{n}=\mathbf{t}\qquad\mathrm{on}\quad\partial\Omega_{N}
\end{equation}
$$

Next, by using the weighted integration and integration by parts, one can obtain:
$$
\begin{equation}
\delta I[\mathbf{u}]=
\int_{\Omega}\sigma_{ij}\delta u_{i,j}dV
-\int_{\Omega}\rho b_{i}\delta u_{i}dV
-\int_{\partial\Omega}t_{i}\delta u_{i}dS=0
\label{eq:weighted-integration}
\tag{2}
\end{equation}
$$
then, the residual for $\eqref{eq:weighted-integration}$ can be expressed as follows:
$$
\begin{equation}
R_{u_{i}}^{I}=
\int_{\Omega}\sigma_{ij}N_{,j}^{I}dV
-\int_{\Omega}\rho b_{i}N^{I}dV
-\int_{\partial\Omega_{N}}t_{i}N^{I}dS
\label{eq:residual}
\tag{3}
\end{equation}
$$
where the "acceleration" term has been ignored.
Then the related jacobian matrix can be given as follows:
$$
\begin{equation}
K_{ik}^{IJ}=\frac{\partial R_{u_{i}}^{I}}{\partial u_{k}^{J}}
=\int_{\Omega}\mathbb{C}\_{ijkl} N_{,j}^{I}N_{,l}^{J}dV
\label{eq:jacobian}
\tag{4}
\end{equation}
$$
where $\mathbb{C}_{ijkl}$ denotes the elasticity tensor. $I$ and $J$ respectively represent the I-th and J-th node of the current element.

## Constitutive laws
For the small strain case, one can have
$$
\begin{equation}
\mathbf{\varepsilon}=\frac{1}{2}(\nabla u+\nabla^{T}u)
\label{eq:strain}
\tag{5}
\end{equation}
$$
which can be easily calculated in AsFem as follows:
```
if(elmtinfo.nDim==1){
  _GradU.SetFromGradU(elmtsoln.gpGradU[1]);
}
else if(elmtinfo.nDim==2){
  _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2]);
}
else if(elmtinfo.nDim==3){
  _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
}
Strain=(_GradU+_GradU.Transpose())*0.5;
```

For the linear elasticity tensor $\mathbb{C}$, you need to give the Young's modulus and poisson ratio, then its value can be defined as follows:
```
Jacobian.SetFromEandNu(E,nu);
```
where AsFem offers the `RankTwoTensor` and `RankFourTensor` class for the complex tensor calculation.

Next, once your strain $\mathbf{\varepsilon}$ and elasticity tensor $\mathbb{C}$ are ready, you can easily get your stress tensor $\mathbf{\sigma}$ as follows:
```
Stress=Jacobian.DoubleDot(Strain);
```

That's all? Yup, all the calculation is done.

If I want to calculate the `vonMises` stress, what should I do? The answer is quite simple, here is the code
```
RankTwoTensor I;
I.SetToIdentity();
devStress=_Stress-_I*(_Stress.Trace()/3.0);
Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
```
once again, the complex tensor calculation is done by our `RankTwoTensor` class.

So, for a general mechanics problem, what you need to tell AsFem are the `stress` and `jacobian` material properties (these two are required by the `mechanics` element) in the following way:
```
Mate.Rank2Materials("stress")=Stress;
Mate.Rank4Materials("jacobian")=Jac;
```
where the values of these two tensors are stored in the Material class.

# Solve the problem

## Define the mesh
We use a square domain here for our calculation, and then the `[mesh]` block can be given as:
```
[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad4
[end]
```
where a $2\times 2$ mesh is defined.

## Gauss point integration (optional)
If one wants to use second order mesh, for instance `meshtype=quad9` or `meshtype=hex27`, then one need to use a higher order gauss points. This can be implemented via:
```
[qpoint]
  type=gauss
  order=3[4]
[end]
```
<span style="color:blue">Normally, you don't need this block!</span>

## Define the DoFs
The DoFs used in this step are the displacements $u_{x}$ and $u_{y}$. Then the `[dofs]` block can be read as:
```
[dofs]
name=ux uy
[end]
```

## Element for the CahnHilliard equation
The model presented in Eq.$\eqref{eq:residual}$ and $\eqref{eq:jacobian}$ can be applied in the following lines
```
[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]
```
where `type=mechanics` option specifies the element we want to use for the linear momentum balance equation. Moreover, we will use the `mymate` material to calculate the necessary material properties required by this element.

## Linear elastic material
By using the following lines in your `[mates]` block, the material properties, i.e., ($\varepsilon$), stress ($\mathbf{\sigma}$), and elasticity tensor ($\mathbb{C}$) can be easily defined:
```
[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3  
    //     E     nu
  [end]
[end]
```
where `type=linearelastic` specifies linear elastic material we want to use. `params=` defines the Youngs modulus $E=210.0$, the poisson ratio $\nu=0.3$.

## Boundary conditions
Now, we want to apply a displacement loading condition to the `top` edge of our rectangle domain and fix the `bottom` edge at the same time, then, we can use:
```
[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux uy
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dofs=uy
    value=0.02
    boundary=top
  [end]
[end]
```


## Static analysis
Again, we need a `[job]` block to start the FEM calculation, which can be given as follows:
```
[job]
  type=static
  debug=dep
[end]
```

## Projection for materials
Hold on for a second, if I want to check the stress, strain, and vonMises stress, what should I do?

The answer is the `[projection]` block, which can do the projection from gauss points to the nodal point.
For example, if we want to save the strain, stress, and vonMises stress into our result file (vtu), one can do:
```
[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]
```
where `scalarmate` and `rank2mate` specify the material name we want to export. Afterwards, they will be displayed in your `Paraview`.




# Run it in AsFem
Now, let's try your first mechanics example in AsFem. You can create a new text file and name it as step6.i or whatever you like. Then copy the following lines into your input file:
```
[mesh]
  type=asfem
  dim=2
  xmax=2.0
  ymax=2.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolids]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3  
    //     E     nu
  [end]
[end]

[bcs]
  [FixUx]
    type=dirichlet
    dofs=ux uy
    boundary=bottom
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dofs=uy
    value=0.02
    boundary=top
  [end]
[end]

[projection]
scalarmate=vonMises
rank2mate=stress strain
[end]

[nonlinearsolver]
  type=nr
  maxiters=20
  r_rel_tol=5.0e-8
  r_abs_tol=4.5e-7
  solver=mumps
[end]

[job]
  type=static
  debug=dep
[end]
```
You can also find the complete input file in `examples/tutorial/step6.i`.



If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):
![ux](ux.jpeg)

![uy](uy.jpeg)

![vonMises](vonMises.jpeg)
