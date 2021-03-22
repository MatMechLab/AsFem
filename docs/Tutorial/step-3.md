---
title: Step 3 - Solve the linear elasticity problem
mathjax: true
date: 2021-01-16 11:09:08
categories:
- Tutorial
tags:
- tutorial
- input file
- mesh
- dofs
- elmts
- mates
- mechanics
- projection
---

# Introduction
In this step, we will try to solve the linear elasticity problem.

# The stress equilibrium equation
The problem we want to solve is the equation of stress equilibrium that reads as follows:
$$
\begin{equation}
\mathbf{\nabla}\cdot\mathbf{\sigma}=\mathbf{0}
\label{eq:stress-eq}
\tag{1}
\end{equation}
$$
where $\mathbf{\sigma}$ denotes Cauchy stress tensor. The constitutive laws for the stress and strain of the small deformation case are below:
$$
\begin{equation}
\mathbf{\sigma}=\mathbb{C}:\mathbf{\epsilon}
\label{eq:stress}
\tag{2}
\end{equation}
$$
and
$$
\begin{equation}
\mathbf{\epsilon}=\frac{1}{2}(\nabla\mathbf{u}+\nabla \mathbf{u}^{T})
\label{eq:dirichlet}
\tag{3}
\end{equation}
$$
where $\mathbf{u}$ is the displacement vector. $\mathbb{C}$ represents the elasticity tensor, which is a function of the Youngs modulus $E$ and Poisson ratio $\nu$.

The related boundary conditions can be read as:
$$
\begin{equation}
\mathbf{t}\cdot\vec{n}=0\qquad\mathrm{on}\quad\partial\Omega_{N}
\label{eq:traction}
\tag{4}
\end{equation}
$$
and
$$
\begin{equation}
\mathbf{u}=\mathbf{u}\_{0}\qquad\mathrm{on}\quad\partial\Omega_{D}
\label{eq:disp}
\tag{5}
\end{equation}
$$
where the traction free condition is assumed in Eq.$\eqref{eq:traction}$

# Define a mesh
For our calculation, we use a rectangular domain here, and then the `[mesh]` block can be given as:
```
[mesh]
  type=asfem
  dim=2
  xmax=5.0
  ymax=5.0
  nx=50
  ny=50
  meshtype=quad4
[end]
```
where a $50\times50$ mesh is defined.


# Define the DoFs
The DoFs used in this step is the displacement vector, namely $u_{x}$ and $u_{y}$. Then the `[dofs]` block can be read as:
```
[dofs]
name=ux uy
[end]
```

# Element for stress equilibrium equation
The model in Eq.$\eqref{eq:stress-eq}$ can be applied in the following lines
```
[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]
```
where `type=mechanics` option specifies the element we want to use for solid mechanics problem. Moreover, we will use the linear elasticity material property, therefore, the related material definition will be given in `mymate` block.

## Linear elasticity material
Via the following lines in your `[mates]` block, the linear elasticity material can be easily defined:
```
[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
  [end]
[end]
```
where `type=linearelastic` specifies linear elasticity material type. `params=` defines the Youngs modulus ($E=210GPa$) and Poisson ratio ($\nu=0.3$)

# Boundary conditions
The boundary conditions, as mentioned in Eq.$\eqref{eq:disp}$, can be applied via the `[bcs]` block. In our case, the *traction* boundary condition in Eq.$\eqref{eq:traction}$ is zero, therefore, only the *displacement* boundary condition needs to be considered:
```
[bcs]
  [fixbottomX]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=bottom
  [end]
  [fixbottomY]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom
  [end]
  [loadY]
    type=dirichlet
    dof=uy
    value=0.1
    boundary=top
  [end]
[end]
```
where we fix the $u_{x}$ and $u_{y}$ to be zero at the bottom edge. While $u_{y}=0.1$ is applied at the top edge of the domain.


# Static analysis
Again, we need a `[job]` block to start the FEM calculation, which can be given as follows:
```
[job]
  type=static
  debug=dep
[end]
```

# Projection
Wait for a minute, where are the stresses and strain? How can I output them? Noooo worries, the `[projection]` block can help you to project the quantities on each `Gauss point` to the `Nodal point` of your mesh. So, if we want to check how the vonMises stress, $\sigma_{xx}$, $\sigma_{yy}$, $\sigma_{xy}$ looks like, we can use:
```
[projection]
name=vonMises stress_xx stress_yy stress_xy
[end]
```
Done!
## Explanation of projection
There is no magic for `[projection]` block, if you take a look at our `MechanicsElmt.cpp` in the `ElmtSystem` class, you will find the following code:
```
gpProj[0]=ScalarMaterials.at("vonMises");
gpProj[1]=Rank2Materials.at("stress")(1,1);//sigma_xx
if(nDim==2){
    gpProj[2]=Rank2Materials.at("stress")(2,2);//sigma_yy
    gpProj[3]=Rank2Materials.at("stress")(1,2);//sigma_xy

    gpProj[4]=Rank2Materials.at("strain")(1,1);//epsilon_xx
    gpProj[5]=Rank2Materials.at("strain")(2,2);//epsilon_yy
    gpProj[6]=Rank2Materials.at("strain")(1,2);//epsilon_xy
}
else if(nDim==3){
    gpProj[2]=Rank2Materials.at("stress")(2,2);//sigma_yy
    gpProj[3]=Rank2Materials.at("stress")(3,3);//sigma_zz
    gpProj[4]=Rank2Materials.at("stress")(2,3);//sigma_yz
    gpProj[5]=Rank2Materials.at("stress")(1,3);//sigma_xz
    gpProj[6]=Rank2Materials.at("stress")(1,2);//sigma_xy

    gpProj[7] =Rank2Materials.at("strain")(1,1);//epsilon_xx
    gpProj[8] =Rank2Materials.at("strain")(2,2);//epsilon_yy
    gpProj[9] =Rank2Materials.at("strain")(3,3);//epsilon_zz
    gpProj[10]=Rank2Materials.at("strain")(2,3);//epsilon_yz
    gpProj[11]=Rank2Materials.at("strain")(1,3);//epsilon_xz
    gpProj[12]=Rank2Materials.at("strain")(1,2);//epsilon_xy
}
```
Yelp, the order of the name you give in your `[projection]` block is related to the `gpProj` quantities in your element. So, if you want to project the quantities in a 3D case, then you will need:
```
[projection]
name=vonMises stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy
[end]
```
Hold on, where is my strain tensor in 2D case? Take it easy, here it is:
```
[projection]
name=vonMises stress_xx stress_yy stress_xy epsilon_xx epsilon_yy epsilon_xy
[end]
```
Done!



# Run it in AsFem
Now, let's try your third example in AsFem. You can create a new text file and name it as step3.i or whatever you like. Then copy the following lines into your input file:
```
[mesh]
  type=asfem
  dim=2
  xmax=5.0
  ymax=5.0
  nx=50
  ny=50
  meshtype=quad4
[end]

[dofs]
name=ux uy
[end]

[elmts]
  [mysolid]
    type=mechanics
    dofs=ux uy
    mate=mymate
  [end]
[end]

[mates]
  [mymate]
    type=linearelastic
    params=210.0 0.3
  [end]
[end]

[bcs]
  [fixbottomX]
    type=dirichlet
    dof=ux
    value=0.0
    boundary=bottom
  [end]
  [fixbottomY]
    type=dirichlet
    dof=uy
    value=0.0
    boundary=bottom
  [end]
  [loadY]
    type=dirichlet
    dof=uy
    value=0.1
    boundary=top
  [end]
[end]

[projection]
name=vonMises stress_xx stress_yy stress_xy
[end]

[job]
  type=static
  debug=dep
[end]
```
You can also find the complete input file in `examples/tutorial/step3.i`.


If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):
![](step3-2d-uy.jpeg)
![](step3-2d-vonMises.jpeg)
