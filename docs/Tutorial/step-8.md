---
title: Step 8 - Write your own model (user-element)
mathjax: true
date: 2022-01-29 20:02:08
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
- umat
- uel
---

# Introduction
In this step, we will try to implement our first user-defined-element (UEL). Now the equation we used here can be whatever you like. Let's take the [heat conduct](https://en.wikipedia.org/wiki/Heat_equation) equation as an example.

The general heat conduct equation takes the following form:
$$
\begin{equation}
\rho c_{p}\frac{\partial T}{\partial t}=\nabla(k\nabla T)+\dot{q}_{v}
\label{eq:T}
\tag{1}
\end{equation}
$$
where $\rho$ is the density, $c_{p}$ denotes the heat capacity, $\dot{q}$ is the volumetric heat source, $k$ represents the thermal conductivity coefficient.

The related boundary conditions can be list below:
$$
\begin{equation}
T=T_{g}~\mathrm{on}~\Omega_{D},\qquad\mathrm{with}~-k\nabla T\cdot\vec{n}=j_{0}~\mathrm{on}~\Omega_{N}
\label{eq:bc}
\tag{2}
\end{equation}
$$

# The user-defined-element (uel)
In AsFem, users can define their own model (governing equations) by using the `uel`. In each uel, user must given the details for the related residual and jacobian calculation.

In this step, we try to write out the code for our heat equation. The residual and system jacobian matrix for Eq.$\eqref{eq:T}$ can be read as follows:

$$
\begin{equation}
\begin{aligned}
R_{T}^{I}&=\int_{\Omega}\rho c_{p}\dot{T} N^{I}dV
+\int_{\Omega}k\nabla T\nabla N^{I}dV
-\int_{\Omega}\dot{q}_{v}N^{I}dV
-\int_{\partial\Omega}k \nabla T\cdot\vec{n}N^{I} dS \\
&=\int_{\Omega}\rho c_{p}\dot{T} N^{I}dV
+\int_{\Omega}k\nabla T\nabla N^{I}dV
-\int_{\Omega}\dot{q}_{v}N^{I}dV
+\int_{\partial\Omega}j_{0}N^{I} dS
\end{aligned}
\label{eq:residual}
\tag{3}
\end{equation}
$$

and
$$
\begin{equation}
K_{TT}^{IJ}=\frac{\partial R_{T}^{I}}{\partial T^{J}}=\int_{\Omega}\rho c_{p}\frac{\partial\dot{T}}{\partial T}N^{J}N^{I}dV
+\int_{\Omega}k\nabla N^{J}\nabla N^{I}dV
\label{eq:jacobian}
\tag{4}
\end{equation}
$$

Then, the formulation part is done!

## Writing code for your uel-1
AsFem offers several uel(1~20), which means you can easily write your code by editing the cpp file in the `src/ElmtSystem` folder. In this case, we will use `uel1`, then one can open the `User1Elmt.cpp` file with whichever text/code editor he/she likes.


### Material  properties
It is not necessary to call a material code for the calculation, but it could be very flexible if one combines `umat` and `uel`. Therefore, one can use one `umat` to calculate the material properties required by Eq.$\eqref{eq:residual}$ and Eq.$\eqref{eq:jacobian}$.

For example, the built-in `thermalmate` defines:
```
Mate.ScalarMaterials("rho")=InputParams[0];// density
Mate.ScalarMaterials("Cp") =InputParams[1];// heat capacity
Mate.ScalarMaterials("K")  =InputParams[2];// thermal conductivity
Mate.ScalarMaterials("Q")  =InputParams[3];// body heat source
```
therefore, you can use `User1Mate` or `UserXMate` for the same purpose. Once your materials are ready, we can move on to the next step, your **first uel**!

### UEL
For Eq.$\eqref{eq:residual}$, namely the residual computation, one can write:
```
localR(1)=_rho*_Cp*soln.gpV[1]*shp.test
         +_K*(soln.gpGradU[1]*shp.grad_test)
         -_Q*shp.test;
```
next, the system jacobian (Eq.$\eqref{eq:jacobian}$) can be calculated as follows:
```
localK(1,1)=_rho*_Cp*shp.trial*shp.test*ctan[1]
           +_K*(shp.grad_trial*shp.grad_test)*ctan[0];
```

That's all the code for your model, done!


# Solve the problem
## Choose the uel-1
Since you have wrote the code for your own element, then you should save the `User1Elmt.cpp` file and **make AsFem again** by running(you should have the **Makefile**, otherwise, please do `cmake CMakeLists.txt`):
```
make -j4
```
Then, you can tell AsFem to use the `User1Elmt` via:
```
[elmts]
  [mythermal]
    type=user1
    dofs=T
    mate=mymate
  [end]
[end]
```
here `type=user1` means we will use the heat equation defined in `User1Elmt.cpp`. The related material properties can be defined as follows:
```
[mates]
  [mymate]
    type=user2
    params=1.0 1.5 2.0 0.0
    //     rho Cp  K   Q
  [end]
[end]
```


## Projection for materials
If one want to check the gradient of the temperature, one can define a vector material as follows:
```
Mate.VectorMaterials("gradT")=elmtsoln.gpGradU[1];
```
then in your `[projection]` block, you can use:
```
[projection]
vectormate=gradT
[end]
```
Then you will see your `gradT` in the `Paraview`.




# Run it in AsFem
Now, let's try your first uel example in AsFem. You can create a new text file and name it as step8.i or whatever you like. Then copy the following lines into your input file:
```
[mesh]
  type=asfem
  dim=2
  xmax=10.0
  ymax=2.0
  nx=100
  ny=20
  meshtype=quad9
[end]

[dofs]
name=T
[end]

[elmts]
  [mythermal]
    type=user1
    dofs=T
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=user2
    params=1.0 1.5 2.0 0.0
    //     rho Cp  K   Q
  [end]
[end]

[bcs]
  [flux]
    type=neumann
    dofs=T
    value=-0.1
    boundary=right
  [end]
[end]

[timestepping]
  type=be
  dt=1.0e-5
  dtmax=1.0e-2
  time=1.0e-1
  optiters=3
  adaptive=true
[end]

[projection]
vectormate=gradT
[end]

[job]
  type=transient
  debug=dep
[end]
```
You can also find the complete input file in `examples/tutorial/step8.i`.
