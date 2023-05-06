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
The problem we want to solve is the equation of stress equilibrium (without body force and acceleration) that reads as follows:
$$
\begin{equation}
\mathbf{\nabla}\cdot\mathbf{\sigma}=\mathbf{0}
\label{eq:stress-eq}
\tag{1}
\end{equation}
$$
where $\mathbf{\sigma}$ denotes Cauchy stress tensor. Below are the constitutive laws for stress and strain in the case of small deformations:
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
We will be using a rectangular domain for our calculations. To define this domain, we use the `mesh` block, which can be expressed as follows:
```
"mesh":{
		"type":"asfem",
		"dim":2,
		"nx":50,
		"ny":50,
		"xmax":5.0,
		"ymax":5.0,
		"meshtype":"quad4",
		"savemesh":true
	},
```
where a $50\times50$ mesh is defined.


# Define the DoFs
In this step, the displacement vector is used as the degree of freedom, specifically $u_{x}$ and $u_{y}$. The `dofs` block can be expressed as follows:
```
"dofs":{
		"names":["ux","uy"]
	}
```

# Element for stress equilibrium equation
The model given in Eq.$\eqref{eq:stress-eq}$ can be implemented using the following lines:
```
"elements":{
		"elmt1":{
			"type":"mechanics",
			"dofs":["ux","uy"],
			"material":{
				"type":"linearelastic",
				"parameters":{
					"E":1.0e3,
					"nu":0.3
				}
			}
		}
	}
```
The `"type":"mechanics"` option in the above lines specifies the element to be used for the solid mechanics problem. Additionally, since we will be using linear elasticity material properties, the related material definition will be provided in the `material` block.

## Linear elasticity material
The following lines within the `material` block can be used to define the linear elasticity material:
```
"material":{
				"type":"linearelastic",
				"parameters":{
					"E":1.0e3,
					"nu":0.3
				}
			}
```
where `"type":"linearelastic"` specifies linear elasticity material model. `parameters` defines the Youngs modulus ($E=100GPa$) and Poisson ratio ($\nu=0.3$). 

##<span style="color:red">**It is worth noting that AsFem does not enforce any specific unit system or interpret results in a particular unit. Therefore, the user must ensure that the input parameters and result interpretation are consistent with their intended application.**</span>.

# Boundary conditions
The boundary conditions mentioned in Eq.$\eqref{eq:disp}$ can be implemented using the `bcs` block. In this case, since the traction boundary condition in Eq.$\eqref{eq:traction}$ is zero, we only need to consider the displacement boundary condition.
```
"bcs":{
		"fix":{
			"type":"dirichlet",
			"dofs":["ux","uy"],
			"bcvalue":0.0,
			"side":["bottom"]
		},
		"load":{
			"type":"dirichlet",
			"dofs":["uy"],
			"bcvalue":0.1,
			"side":["top"]
		}
	}
```
where we fix $u_{x}$ and $u_{y}$ to be zero at the bottom edge of the domain, while $u_{y}=0.1$ is applied at the top edge.


# Static analysis
To start the FEM calculation, we need to use the `job` block. The following lines can be used for this purpose:
```
"job":{
		"type":"static",
		"print":"dep"
	}
```

# Projection
Wait for a second, where are the stresses and strain? How can I output them? 

Noooo worries, the `projection` block can be used to project the quantities from each `Gauss point` to the `Nodal point` of your mesh. To visualize the vonMises stress, $\sigma_{xx}$, $\sigma_{yy}$, and $\sigma_{xy}$, you can use the following code:
```
"projection":{
		"type":"default",
		"scalarmate":["vonMises-stress"],
		"rank2mate":["stress","strain"]
	}
```
The default option for the `projection` block (`"type":"default"`) will execute the simplified least squares projection algorithm.


## Explanation of projection
The `projection` block is not a magic tool. If you examine the LinearElasticMaterial.cpp file in the MateSystem class, you will find the following code:
```
mate.ScalarMaterial("vonMises-stress")=sqrt(1.5*m_devStress.doubledot(m_devStress));
mate.Rank2Material("strain")=m_strain;
mate.Rank2Material("stress")=m_stress;
```
The material properties, which include scalar, vector, and tensor materials, are stored in the material class and can be easily accessed using their respective material names. Users only need to provide the names of the output material properties in the `projection` block. Additionally, the general tensor form used in this block is designed to work for all dimensions, from 1D to 3D.



# Run it in AsFem
To try the third example in AsFem, create a new text file and name it `step3.json` or any other name you prefer. Then, copy and paste the following lines into your input file:
```
{
	"mesh":{
		"type":"asfem",
		"dim":2,
		"nx":50,
		"ny":50,
		"xmax":5.0,
		"ymax":5.0,
		"meshtype":"quad4",
		"savemesh":true
	},
	"dofs":{
		"names":["ux","uy"]
	},
	"elements":{
		"elmt1":{
			"type":"mechanics",
			"dofs":["ux","uy"],
			"material":{
				"type":"linearelastic",
				"parameters":{
					"E":1.0e3,
					"nu":0.3
				}
			}
		}
	},
	"projection":{
		"type":"default",
		"scalarmate":["vonMises-stress"],
		"rank2mate":["stress","strain"]
	},
	"nlsolver":{
		"type":"newton",
		"solver":"gmres",
		"maxiters":50,
		"abs-tolerance":5.0e-7,
		"rel-tolerance":5.0e-10,
		"s-tolerance":0.0,
		"preconditioner":"lu"
	},
	"output":{
		"type":"vtu",
		"interval":1
	},
	"bcs":{
		"fix":{
			"type":"dirichlet",
			"dofs":["ux","uy"],
			"bcvalue":0.0,
			"side":["bottom"]
		},
		"load":{
			"type":"dirichlet",
			"dofs":["uy"],
			"bcvalue":0.1,
			"side":["top"]
		}
	},
	"qpoints":{
		"bulk":{
			"type":"gauss-legendre",
			"order":2
		}
	},
	"job":{
		"type":"static",
		"print":"dep"
	}
}
```

<span style="color:red">It's worth noting that if the `"preconditioner":"lu"` option in the `nlsolver` block doesn't work, you may need to either change `lu` to `jacobi`, or install MUMPS or SuperLU for your PETSc</span>.

You can also find the complete input file in `examples/tutorial/step3-2d.json`.


If everything goes well, you can see the following images in your [Paraview](https://www.paraview.org/download/):
![](step3-2d-uy.jpeg)
![](step3-2d-vonMises.jpeg)
