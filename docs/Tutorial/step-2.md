

# Introduction
In the first step of our analysis, which is detailed in [step-1](https://m3group.github.io/AsFem/Tutorial/step-1/), we defined the mesh. However, additional information is required to conduct a basic finite element analysis. As such, we will proceed by defining the degrees of freedom (DoFs) and the model in this step. Finally, we will apply the relevant boundary conditions to obtain the desired solution.

# The poisson equation
The issue we want to solve is the linear Poisson equation that reads as follows:
$$
\begin{equation}
\sigma\nabla^{2}\phi=f
\label{eq:poisson}
\tag{1}
\end{equation}
$$
where $\sigma$ and $f$ denote the model's coefficients. Below are the related boundary conditions:
$$
\begin{equation}
\sigma\nabla\phi\cdot\vec{n}=0\qquad\mathrm{on}\quad\partial\Omega_{N}
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
where the *Dirichlet* boundary condition and *Neumann* boundary condition are represented by the subscripts $D$ and $N$.


# Define the degree of freedom (DoF)

The degree of freedom (DoF) or the degrees of freedom (DoFs) can be used to define the name of each DoF and also to apply the necessary boundary conditions (`[bcs]`), elements (`[elmts]`), and so on. The `[dofs]` block looks like below:

The degrees of freedom (DoF), or DoFs, are used to assign names to each DoF and to apply the required boundary conditions (`bcs`) and elements (`elements`). The block of code defining the DoFs is presented below as `dofs`:
```
"dofs":{
		"names":["dof-a","dof-b","dof-c",...,"dof-x"]
	}
```
## Options
The `dofs` block includes an option called names, which is used to specify the name of each degree of freedom. It is important to note that the order in which the names are listed corresponds to the index of each degree of freedom. For example, in a 2D elastic analysis, two displacement degrees of freedom are required, namely `disp_x` and `disp_y`. The `dofs` block should thus be defined as follows:
```
"dofs":{
		"names":["disp_x","disp_y"]
	}
```
In this specific case, `disp_x` is assigned as the first degree of freedom with an index of 1, while `disp_y` is assigned as the second degree of freedom with an index of 2. Note that the only option available in the `dofs` block is `names`, with no additional options

Since the Poisson equation involves only a single degree of freedom, the final expression for the `dofs` block is as follows
```
"dofs":{
		"names":["phi"]
	}
```

# Element for Poisson equation
With the degree of freedom established, we now need to define the model in Eq.$\eqref{eq:poisson}$. To accomplish this, we use the `elements` block, which is defined as follows:
```
"elements":{
		"elmt1":{
			"type":"poisson",
			"dofs":["phi"],
			"material":{
				"type":"constpoisson",
				"parameters":{
					"sigma":1.0,
					"f":1.0e1
				}
			}
		}
	}
```
The `type` option specifies the element or model to be used, which may be one of the built-in elements in AsFem or a user-defined element (`UEL`). The degrees of freedom to be used in the element are defined by the `dofs` option, while the material option specifies the name of the material sub-block to be used. Once the `elements` block is defined, the model described in Eq.$\eqref{eq:poisson}$ is complete.

## Material properties
The coefficients $k$ and $F$, which correspond to the material properties in Eq.$\eqref{eq:poisson}$, can be calculated or defined using the `material` sub-block, as shown below:
```
"material":{
				"type":"constpoisson",
				"parameters":{
					"sigma":1.0,
					"f":1.0e1
				}
			}
```
The `material` sub-block includes the `type` option, which specifies the name of the material type defined in AsFem. The parameters option is used to define the parameters that will be used in our model. In this case, we will use $\sigma=1.0$ and $f=10.0$

# Boundary conditions
The boundary conditions described in Eq.$\eqref{eq:dirichlet}$ and $\eqref{eq:neumann}$ can be applied using the `bcs` block. In our case, the Neumann boundary condition in Eq.$\eqref{eq:neumann}$ is zero, and thus, only the Dirichlet boundary condition needs to be taken into account:
```
"bcs":{
		"left":{
			"type":"dirichlet",
			"dofs":["phi"],
			"bcvalue":0.1,
			"side":["left"]
		},
		"right":{
			"type":"dirichlet",
			"dofs":["phi"],
			"bcvalue":0.5,
			"side":["right"]
		}
	}
```
The `left` and `right` sub-block is used to specify different boundary conditions for the left and right boundaries. The `type` option within each sub-block is used to specify the type of boundary condition, which is supported by AsFem. The `dof` option denotes the name of the DoF to which the given boundary conditions will be applied. In our case, we constrain the value of $\phi$ on the left and right sides of a rectangular domain to be 0.1 and 0.5, respectively.


# Static analysis
Now that the model and boundary conditions have been defined, we are ready to begin the FEM calculation. To specify the type of analysis we want to perform, we use the `job` block. In this case, since we are conducting a static analysis, the block can be defined as follows:
```
"job":{
		"type":"static",
		"print":"on",
	}
```
If you want to observe how the iteration information changes during the analysis, you can use the following:
```
"job":{
		"type":"static",
		"print":"dep",
		"restart":true
	}
```
The `print` option is used to enable basic information output in the terminal. If you prefer to see fewer outputs, you can use `"print":"off"`.

# Run it in AsFem
To proceed with the second example in AsFem, create a new text file and name it `step2.json`, or a name of your choosing. Next, copy the following lines into your input file:
```
{
	"mesh":{
		"type":"asfem",
		"dim":2,
		"nx":50,
		"ny":50,
		"meshtype":"quad4",
		"savemesh":true
	},
	"dofs":{
		"names":["phi"]
	},
	"elements":{
		"elmt1":{
			"type":"poisson",
			"dofs":["phi"],
			"material":{
				"type":"constpoisson",
				"parameters":{
					"sigma":1.0,
					"f":1.0e1
				}
			}
		}
	},
	"projection":{
		"type":"default",
		"vectormate":["gradu"]
	},
	"bcs":{
		"left":{
			"type":"dirichlet",
			"dofs":["phi"],
			"bcvalue":0.1,
			"side":["left"]
		},
		"right":{
			"type":"dirichlet",
			"dofs":["phi"],
			"bcvalue":0.5,
			"side":["right"]
		}
	},
	"nlsolver":{
		"type":"newton",
		"solver":"gmres",
		"maxiters":50,
		"abs-tolerance":5.0e-7,
		"rel-tolerance":5.0e-10,
		"s-tolerance":0.0,
		"preconditioner":"jacobi"
	},
	"output":{
		"type":"vtu",
		"interval":1
	},
	"qpoints":{
		"bulk":{
			"type":"gauss-legendre",
			"order":2
		}
	},
	"job":{
		"type":"static",
		"print":"dep",
		"restart":true
	}
}
```
You can also find the complete input file in `examples/tutorial`.


If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):
![](step2.jpeg)



Wait a second, what should I do if I want to solve a Poisson equation in 3D?


The solution is straightforward - you just need to modify your mesh to be 3D. Here is an example of how you can do this using the input file `step2-3d.json`:
```
"mesh":{
		"type":"asfem",
		"dim":3,
		"nx":50,
		"ny":50,
		"nz":50,
		"meshtype":"hex8",
		"savemesh":true
	}
```
then you will see:
![](step2-3d.jpeg)
