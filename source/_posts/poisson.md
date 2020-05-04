---
title: Poisson equation
mathjax: true
date: 2020-05-03 16:08:45
categories:
- Examples
tags:
- Examples
- Poisson
---

A classical Poisson equation can be written out as follow:
$$\begin{equation}
\nabla\cdot(\sigma(\phi)\nabla\phi)=F(\phi)
\end{equation}$$
with the boundary conditions as follows:
$$\begin{align}
\sigma(\phi)\nabla\phi\cdot\vec{n}&=0\qquad\mathrm{on}\quad\partial\Omega_{N}\\
\phi&=\phi_{0}\qquad\mathrm{on}\quad\partial\Omega_{D}
\end{align}$$
where both the $\sigma(\phi)$ and $F(\phi)$ can be the function of $\phi$. In this scenario, the poisson equation becomes the nonlinear case.

By choosing a suitable test function $\delta\phi$, we can have the integration of the equation as follow:
$$\begin{equation}
\int_{\Omega}\nabla\cdot(\sigma(\phi)\nabla\phi)\delta\phi dV=\int_{\Omega}F(\phi)\delta\phi dV
\end{equation}$$
applying the divergence theorem and integrating by parts, one can have the weak form as follow:
$$\begin{equation}
\int_{\partial\Omega}\sigma(\phi)\nabla\phi\cdot\vec{n}\delta\phi dS
-\int_{\Omega}\sigma(\phi)\nabla\phi\nabla\delta\phi dV
-\int_{\Omega}F(\phi)\delta\phi dV=0
\end{equation}$$

Then, by applying the boundary conditions, the residual can be stated as:
$$\begin{equation}R_{\phi}^{I}=
\int_{\Omega}\sigma(\phi)\nabla\phi\nabla N^{I}dV
+\int_{\Omega}F(\phi)N^{I}dV
\end{equation}$$
where the superscript $I$ denotes the i-th node of the current element, and where the subscript $\phi$ represents the related degree of freedom (DoF).

Furthermore, one can have the tangential matrix or the stiffness matrix as follows:
$$\begin{equation}
\begin{aligned}
K_{\phi\phi}^{IJ}=\frac{\partial R^{I}_{\phi}}{\partial\phi^{J}}
&=\int_{\Omega}\frac{\partial\sigma(\phi)}{\partial\phi}N^{J}\nabla\phi\nabla N^{I}dV
+\int_{\Omega}\sigma(\phi)\nabla N^{J}\nabla N^{I}dV\\
&+\int_{\Omega}\frac{\partial F(\phi)}{\partial\phi}N^{J}N^{I}dV
\end{aligned}
\end{equation}
$$


Remember, in order to make our code as simple as possible, the derivative of K is without the "-" sign(according to the mathematical definition of Newton-Raphson, we should obtain the K as: $K=-\partial R/\partial\phi$). However, in AsFem, we ignore the "-" sign to make the programming as easy as possible for our users.

Then the one-to-one mapping between the code and the formula can be written out as follow:

```c++
if(isw==3||isw==6){
    for(int i=1;i<=nNodes;++i){
        rhs(i)+=ScalarMaterials[0]*(gpGradU[0]*shp.shape_grad(i))
                   +ScalarMaterials[2]*shp.shape_value(i);
        if(isw==6){
            for(int j=1;j<=nNodes;++j){
                K(i,j)+=ScalarMaterials[1]*shp.shape_value(j)*(gpGradU[0]*shp.shape_grad(i))*ctan[0]
                       +ScalarMaterials[0]*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]
                       +ScalarMaterials[3]*shp.shape_value(j)*shp.shape_value(i)*ctan[0];
            }
        }
    }
}
```

the ScalarMaterials used here is calculated from the Materials in AsFem, for the details one is referred to [ConstPoissonMate](https://github.com/yangbai90/AsFem/blob/master/src/MateSystem/ConstPoissonMaterial.cpp)
```c++
// In constpoisson material:
// _MateValues[0]=1.0;// sigma
// _MateValues[1]=0.0;// dsigma/dphi
// _MateValues[2]=1.0;// F
// _MateValues[3]=0.0;// dF/dphi
```


In order to solve this equation, we should define a domain for it:
```
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  nx=100
  ny=100
  meshtype=quad4
[end]
```
where the $100\times100$ mesh for a $10\times10$ domain is applied.

Then we need to define the name of our Dof:
```
[dofs]
name=phi
[end]
```
for sure, the name can be whatever you like.

Next, we need to tell AsFem we want to use the element that can do the calculation for the equations we mentioned above:
```
[elmts]
    [poisson]
	  type=poisson
	  dofs=phi
	  mate=linear
      domain=alldomain
    [end]
[end]
```
Here the 'type=poisson' tells AsFem that we will use the element for the Poisson equation, and it will use the information of 'dof=phi' for the calculation. Later, the material we used in this element will look for the material block '[linear]' since we tell him the 'mate=linear'. Finally, the equation is applied to the domain 'alldomain'(which is the default value for the whole domain, and it can be ignored.)

Then, after the [elmts] block is defined, we can tell AsFem, how the 'linear' material block should look like:
```
[mates]
  [linear]
    type=constpoisson
    params=1.0e1 0.3
  [end]
[end]
```
here the 'type=constpoisson' will ask AsFem to run the constpoisson material, and the parameters we used for $\sigma$ and $F$ are 1.0e1 and 0.3.

Now, all the equations we need are ready to get the final solution, we need to apply the correct boundary conditions:
```
[bcs]
  [fixphi]
    type=dirichlet
    dof=phi
    boundary=left right bottom top
    value=0.0
  [end]
[end]
```
Here we fixed the four edges of the plane to be 0.0. If you only want the top and bottom edge to be fixed, then you just need to give: 'boundary=top bottom'.

Now we can assign a job to AsFem and ask him to give us the final result:
```
[job]
  type=static
  debug=dep
[end]
```
where 'debug=dep' or 'debug=true' allows AsFem to print out some necessary information. For sure, you can also disable it via: 'debug=false'.

If we run it, the result should look like:
![](poisson.jpeg)

The complete input file can be found here: [input](https://github.com/yangbai90/AsFem/blob/master/tests/poisson/quad4_linear.i)
