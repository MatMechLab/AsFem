---
title: Phasefield fracture model in chemo-mech coupled case
author: Mengmeng Li
mathjax: true
date: 2025-07-08 16:02:08
categories:
- Examples
tags:
- examples
- input file
- mesh
- dofs
- elmts
- mates
- stress
- umat
- fracture
---

---
<span style="color: #6B7280; font-weight: bold;">Author: Mengmeng Li</span> • [MMlab 🔗](https://www.x-mol.com/groups/matmechlab)

# Introduction
In this example, we will try to solve the coupling diffusion and Allen–Cahn Fracture model.


# Model
In this model, the damage phase is described by an order parameter $d$ which varies smoothly from 0 (undamaged case) to 1 (fully damaged case). Therefore, the system free energy can be read as follows:

The fracture model follows Allen–Cahn Fracture equation:
$$
\begin{equation}
\frac{\partial \eta}{\partial t} = -L \frac{\delta \psi}{\delta \eta}
\label{eq:AC_equation}
\tag{1}
\end{equation}
$$

The diffusion model is list as follows:
$$
\begin{equation}
\frac{\partial c}{\partial t} = \nabla \cdot (D \nabla c -D c\Omega \nabla{\sigma_{h}})
\label{eq:diffusion_equation}
\tag{2}
\end{equation}
$$
where $\sigma_{h}$ is hrdyostatic stress.
 
The governing equaitons for this model are list below:
$$
\begin{equation}
\mathbf{\nabla}\cdot\boldsymbol{\sigma}=\boldsymbol{0}
\label{eq:stress-equilibrium}
\tag{3}
\end{equation}
$$
<!-- The driving force for damage development is listed below: -->
<!--
$$
\begin{equation}
\frac{\partial d}{\partial t}=2(1-d)\mathcal{H}-\frac{\mathcal{G}_{c}}{l}(d-l^{2}\Delta d)
\label{eq:damage-equation}
\tag{4}
\end{equation}
$$ -->

# Input file
In the example/pffracture folder, there are several `geo` file, which can be used to generate the `msh` mesh file for your simulation. Before we start, you should have the mesh file.

## Mesh
In this case, we'll import the mesh from `gmsh`, which can be done as follows:
```
"mesh":{
  "type":"msh2",
  "file":"rect50grs.msh",
  "savemesh":true
},
```
Here, you should use [gmsh](https://gmsh.info/) to generate the related msh file!

## Dofs
Next, you need to define the dofs, it should be $d$, $u_{x}$, and $u_{y}$ for 2d case, and $d$, $u_{x}$, $u_{y}$, and $u_{z}$ for 3d case. The `[dofs]` block should looks like:
```
"dofs":{
  "names":["c","d","ux","uy"]
},
```
```
"dofs":{
  "names":["c","d","ux","uy","uz"]
},
```
Here, the sequence of your `dofs` name matters !!!

## elmts and mates
The governings in Eq.$\eqref{eq:stress-equilibrium}$ and Eq.$\eqref{eq:damage-equation}$ can be implemented by using the following `elmts`:
```
"elements":{
  "elmt1":{
    "type":"diffusionacfracture",
    "dofs":["c","d","ux","uy"],
    "domain":["1"],
    "material":{
      "type":"diffusionacfracture",
      "parameters":{
        "stabilizer":1.0e-5,
        "L":1.0e4,
        "Gc":2.7e-3,
        "eps":0.01,
        "K":121.15,
        "G":80.77,
        "Cref":0.0,
        "D":0.1,
        "Omega": 5.65895e-02
      }
    }
  },
```
The parameter `type=diffusionacfracture` instructs AsFem to call the diffusion and Allen-Cahn fracture model. This example demonstrates a polycrystalline model where `elmt1` represents the first grain, , with each grain assigned a different value of $\Omega$. The complete input file is available at `examples/pffracture/rect50grsdifffrac.json`.

As mentioned previously, several material parameters are essential for the calculations: the bulk modulus ($K$), shear modulus ($G$), fracture energy ($\mathcal{G}_{c}$), and the regularization length ($l$, also denoted as `eps`) which governs damage evolution.

## boundary condition
For the tensile test, one can use the following `bcs`:
```
  "bcs":{
    "fixux":{
      "type":"dirichlet",
      "dofs":["ux"],
      "bcvalue":0.0,
      "side":["left"]
    },
    "fixuy":{
      "type":"dirichlet",
      "dofs":["uy"],
      "bcvalue":0.0,
      "side":["bottom"]
    },
    "load":{
      "type":"neumann",
      "dofs":["c"],
      "bcvalue":"-0.015",
      "side":["top","right"]
    }
  },
```

# Run it in AsFem
Now, let's try your first finite strain fracture model in AsFem. You can create a new text file and name it as `diffusionacfracture.json` or whatever you like. A part of input file is as follows:
```
{
  "mesh":{
    "type":"msh2",
    "file":"rect50grs.msh",
    "savemesh":true
  },
  "dofs":{
    "names":["c","d","ux","uy"]
  },
  "elements":{
    "elmt1":{
      "type":"diffusionacfracture",
      "dofs":["c","d","ux","uy"],
      "domain":["1"],
      "material":{
        "type":"diffusionacfracture",
        "parameters":{
          "stabilizer":1.0e-5,
          "L":1.0e4,
          "Gc":2.7e-3,
          "eps":0.01,
          "K":121.15,
          "G":80.77,
          "Cref":0.0,
          "D":0.1,
          "Omega":   5.65895e-02
        }
      }
    },
  
  ...............................................................
  ...............................................................
  },

  "projection":{
    "type":"default",
    "scalarmate":["vonMises-stress","vonMises-strain","hydrostatic-stress"],
    "rank2mate":["stress"]
  },
  "linearsolver":{
    "type":"mumps",
    "preconditioner":"lu",
    "maxiters":10000,
    "restarts":2600,
    "tolerance":1.0e-26
  },
  "nlsolver":{
    "type":"newton",
    "maxiters":100,
    "abs-tolerance":7.5e-7,
    "rel-tolerance":1.0e-10,
    "s-tolerance":0.0
  },
  "output":{
    "type":"vtu",
    "interval":10
  },
  "bcs":{
    "fixux":{
      "type":"dirichlet",
      "dofs":["ux"],
      "bcvalue":0.0,
      "side":["left"]
    },
    "fixuy":{
      "type":"dirichlet",
      "dofs":["uy"],
      "bcvalue":0.0,
      "side":["bottom"]
    },
    "load":{
      "type":"neumann",
      "dofs":["c"],
      "bcvalue":"-0.015",
      "side":["top","right"]
    }
  },
  "qpoints":{
    "bulk":{
      "type":"gauss-legendre",
      "order":2
    }
  },
  "timestepping":{
    "type":"be",
    "dt0":1.0e-6,
    "dtmax":1.0e-1,
    "dtmin":1.0e-12,
    "optimize-iters":5,
    "end-time":4.8e1,
    "growth-factor":1.1,
    "cutback-factor":0.85,
    "adaptive":true
  },
  "postprocess":{
    "vonMises-stress":{
      "type":"nodalscalarmate",
      "parameters":{
        "nodeid":662,
        "scalarmate":"vonMises-stress"
      }
    },
    "uy":{
      "type":"nodalvalue",
      "dof":"uy",
      "parameters":{
        "nodeid":662
      }
    }
  },
  "job":{
    "type":"transient",
    "print":"dep"
  }
}

```
Due to space constraints, only element1 is displayed here. For the complete input file, please refer to `examples/pffracture/rect50grsdifffrac.json`.


If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):

![diffusionacFracture](diffusionacFracture.png)