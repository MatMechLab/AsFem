---
title: Phasefield fracture model in finite strain case
author: Mengmeng Li
date: 2025-07-08
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
- phase field
- fracture
---
---

<!-- **Author: Mengmeng Li** • [MMlab 🔗](https://www.x-mol.com/groups/matmechlab) -->
<span style="color: #6B7280; font-weight: bold;">Author: Mengmeng Li</span> • [MMlab 🔗](https://www.x-mol.com/groups/matmechlab)

# Introduction
In this example, we will try to solve the phase-field finite strain fracture model.

# Model
In this model, the damage phase is described by an order parameter $d$ which varies smoothly from 0 (undamaged case) to 1 (fully damaged case). Therefore, the system free energy can be read as follows:
$$
\begin{equation}
\psi=\psi_{\mathrm{ela}}+\psi_{\mathrm{frac}}
\label{eq:psi}
\tag{1}
\end{equation}
$$
where
$$
\begin{equation}
\psi_{\mathrm{ela}} = g(d)\psi_{\mathrm{ela}}^{+}+\psi_{\mathrm{ela}}^{-}
\tag{2}
\end{equation}
$$
where $g(d)$ is the degradation function, and the positive energy $\psi_{\mathrm{ela}}^{+}$ and the negative energy $\psi_{\mathrm{ela}}^{-}$ are given as follows:

$$
\begin{equation}
\psi_{\mathrm{ela}}^{+}=
\begin{cases}
\frac{K}{2}(\frac{J_e^2-1}{2}-\ln(J_e))+\frac{\mu}{2}(\bar{I}_1-3) & \mathrm{if} J_e \geq 1 \\
\frac{\mu}{2}(\bar{I}_1-3) & \mathrm{if} J_e < 1
\end{cases}
\tag{3}
\label{eq:psi-e-positive}
\end{equation}
$$

and

$$
\begin{equation}
\psi_{\mathrm{ela}}^{-}=
\begin{cases}
0 & \mathrm{if } J_e \geq 1 \\
\frac{K}{2}\left(\frac{J_e^2-1}{2}-\ln(J_e)\right) & \mathrm{if } J_e < 1
\end{cases}
\tag{4}
\label{eq:psi-e-negative}
\end{equation}
$$

Thus the positive and negative part of the 2nd PK stress $\boldsymbol{S}$ can take the following
forms:

$$
\begin{equation}
\boldsymbol{S}^+=
\frac{\partial \psi_{\text{ela}}^+}{\partial \boldsymbol{E}_e} = 2 \frac{\partial \psi_{\text{ela}}^+}{\partial\boldsymbol{C}_e} = \begin{cases}
\frac{K}{2}(J_e^2-1)\boldsymbol{C}_e^{-1}+\mu J_e^{-\frac{2}{3}}(I-\frac{1}{3}I_1\boldsymbol{C}_e^{-1}) & \mathrm{if } J_e \geq 1 \\
\mu J_e^{-\frac{2}{3}}(\boldsymbol{I} - \frac{1}{3}I_1 \boldsymbol{C}_e^{-1}) & \text{if } J_e < 1
\end{cases}
\tag{5}
\label{eq:stress-e-positive}
\end{equation}
$$

and 

$$
\begin{equation}
\boldsymbol{S}^-=
\frac{\partial \psi_{\text{ela}}^-}{\partial \boldsymbol{E}_e} = 2 \frac{\partial \psi_{\text{ela}}^-}{\partial \boldsymbol{C}_e} = \begin{cases}
0 & \text{if } J_e \geq 1 \\
\frac{K}{2}(J_e^2 - 1)\boldsymbol{C}_e^{-1} & \text{if } J_e < 1
\end{cases}
\tag{6}
\label{eq:stress-e-negative}
\end{equation}
$$

The fracture free energy can be given by:

$$
\begin{equation}
\psi_{\mathrm{frac}}=\mathcal{G}_{c}(\frac{d^{2}}{2l}+\frac{l}{2}\lvert\nabla d\rvert^{2})
\label{eq:psi-frac}
\tag{7}
\end{equation}
$$
with $\mathcal{G}_{c}$ and $l$ being the critical energy release rate and the length scale parameter, respectively. 

The governing equaitons for this model are list below:
$$
\begin{equation}
\mathbf{\nabla}\cdot\boldsymbol{P}=\boldsymbol{0}
\label{eq:stress-equilibrium}
\tag{8}
\end{equation}
$$
and
$$
\begin{equation}
\frac{\partial d}{\partial t}=2(1-d)\mathcal{H}-\frac{\mathcal{G}_{c}}{l}(d-l^{2}\Delta d)
\label{eq:damage-equation}
\tag{9}
\end{equation}
$$
The history variable $\mathcal{H}$ is calculated as follows:

$$
\begin{equation}
\mathcal{H}=
\begin{cases}
\psi_{e}^{+} & \mathrm{if}\quad\psi_{e}^{+}>\mathcal{H}_{n}\\
\mathcal{H}_{n} &\mathrm{otherwise}
\end{cases}
\label{eq:hist}
\tag{10}
\end{equation}
$$

# Input file
In the example/pffracture folder, there are several `geo` file, which can be used to generate the `msh` mesh file for your simulation. Before we start, you should have the mesh file.

## Mesh
In this case, we'll import the mesh from `gmsh`, which can be done as follows:
```
"mesh":{
    "type":"msh2",
    "file":"Tensile3D.msh",
    "savemesh":true
  },
```
Here, you should use [gmsh](https://gmsh.info/) to generate the related msh file!

## Dofs
Next, you need to define the dofs, it should be $d$, $u_{x}$, and $u_{y}$ for 2d case, and $d$, $u_{x}$, $u_{y}$, and $u_{z}$ for 3d case. The `[dofs]` block should looks like:
```
"dofs":{
  "names":["d","ux","uy","uz"]
},
```
Here, the sequence of your `dofs` name matters !!!

## elmts and mates
The governings in Eq.$\eqref{eq:stress-equilibrium}$ and Eq.$\eqref{eq:damage-equation}$ can be implemented by using the following `elmts`:
```
"elements":{
  "elmt1":{
    "type":"allencahnfracture",
    "dofs":["d","ux","uy","uz"],
    "material":{
      "type":"neohookeanpffracture",
      "parameters":{
        "L":1.0e6,
        "Gc":2.7e-3,
        "eps":0.012,
        "K":121.15,
        "G":80.77,
        "stabilizer":1.0e-6,
        "finite-strain":true
      }
    }
  }
},
```
where `type=neohookeanpffracture` tells AsFem, our users want to call the neohookean fracture model.

For the material calculation, as mentioned before, several parameters are required, i.e., $E$, $\nu$ for the Youngs modulus and poisson ration, $\eta$, $\mathcal{G}_{c}$, and $l$ for the damage evolution. Therefore, the related `material` block can be given as below.

then, your $\mathcal{H}$ will always use the $\mathcal{H}_{old}$ value from previous step.

## boundary condition
For the tensile test, one can use the following `bcs`:
```
"bcs":{
  "fixux":{
    "type":"dirichlet",
    "dofs":["ux"],
    "bcvalue":0.0,
    "side":["bottom","top"]
  },
  "fixuy":{
    "type":"dirichlet",
    "dofs":["uy"],
    "bcvalue":0.0,
    "side":["bottom"]
  },
  "fixuz":{
    "type":"dirichlet",
    "dofs":["uz"],
    "bcvalue":0.0,
    "side":["bottom","top"]
  },
  "loading":{
    "type":"dirichlet",
    "dofs":["uy"],
    "bcvalue":"0.1*t",
    "side":["top"]
  }
},
```
Done? Yep, all the things is done!


# Run it in AsFem
Now, let's try your first finite strain fracture model in AsFem. You can create a new text file and name it as newhookeanFracture.json or whatever you like. Then copy the following lines into your input file:
```
{
  "mesh":{
    "type":"msh2",
    "file":"Tensile3D.msh",
    "savemesh":true
  },
  "dofs":{
    "names":["d","ux","uy","uz"]
  },
  "elements":{
    "elmt1":{
      "type":"allencahnfracture",
      "dofs":["d","ux","uy","uz"],
      "material":{
        "type":"neohookeanpffracture",
        "parameters":{
          "L":1.0e6,
          "Gc":2.7e-3,
          "eps":0.012,
          "K":121.15,
          "G":80.77,
          "stabilizer":1.0e-6,
          "finite-strain":true
        }
      }
    }
  },
  "projection":{
    "type":"default",
    "scalarmate":["vonMises-stress"],
    "rank2mate":["Stress"]
  },
  "bcs":{
    "fixux":{
      "type":"dirichlet",
      "dofs":["ux"],
      "bcvalue":0.0,
      "side":["bottom","top"]
    },
    "fixuy":{
      "type":"dirichlet",
      "dofs":["uy"],
      "bcvalue":0.0,
      "side":["bottom"]
    },
    "fixuz":{
      "type":"dirichlet",
      "dofs":["uz"],
      "bcvalue":0.0,
      "side":["bottom","top"]
    },
    "loading":{
      "type":"dirichlet",
      "dofs":["uy"],
      "bcvalue":"0.1*t",
      "side":["top"]
    }
  },
  "linearsolver":{
    "type":"mumps",
    "preconditioner":"lu",
    "maxiters":10000,
    "restarts":3600,
    "tolerance":1.0e-26
  },
  "nlsolver":{
    "type":"newton",
    "maxiters":35,
    "abs-tolerance":8.5e-7,
    "rel-tolerance":5.0e-10,
    "s-tolerance":0.0
  },
  "timestepping":{
    "type":"be",
    "dt0":1.0e-5,
    "dtmax":5.0e-3,
    "dtmin":1.0e-12,
    "optimize-iters":4,
    "end-time":5.0e1,
    "growth-factor":1.05,
    "cutback-factor":0.85,
    "adaptive":true
  },
  "output":{
    "type":"vtu",
    "interval":10
  },
  "qpoints":{
    "bulk":{
      "type":"gauss-legendre",
      "order":2
    }
  },
  "postprocess":{
    "uy":{
      "type":"sideaveragevalue",
      "dof":"uy",
      "side":["top"]
    },
    "fy":{
      "type":"sideintegralrank2mate",
      "side":["top"],
      "parameters":{
        "rank2mate":"Stress",
        "i-index":2,
        "j-index":2
      }
    }
  },
  "job":{
    "type":"transient",
    "print":"dep"
  }
}

```
You can also find the complete input file in `examples/pffracture/neohookeanfrac-3d-tensile-check.json`.



If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):

![newHookeanFracture](newHookeanFracture.png)
