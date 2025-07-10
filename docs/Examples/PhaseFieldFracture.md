---
title: Phasefield fracture model in small strain case
mathjax: true
date: 2021-12-24 16:02:08
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
<span style="color: #6B7280; font-weight: bold;">Author: Yang Bai</span> • [MMlab 🔗](https://www.x-mol.com/groups/matmechlab)

# Introduction
In this example, we will try to solve the phase-field fracture model which is implemented based on Prof. Miehe's [Model](https://doi.org/10.1002/nme.2861).


# Model
In this model, the damage phase is described by an order parameter $d$ which varies smoothly from 0 (undamaged case) to 1 (fully damaged case). Therefore, the system free energy can be read as follows:
$$
\begin{equation}
\psi=\psi_{d}+\psi_{e}
\label{eq:psi}
\tag{1}
\end{equation}
$$
where $g(d)$ is the degradation function, and the damage phase free energy $\psi_{d}$ is given as follows:
$$
\begin{equation}
\psi_{d}=\mathcal{G}_{c}(\frac{d^{2}}{2l}+\frac{l}{2}\lvert\nabla d\rvert^{2})
\label{eq:psi-d}
\tag{2}
\end{equation}
$$
with $\mathcal{G}_{c}$ and $l$ being the critical energy release rate and the length scale parameter, respectively.

The elastic free energy can be defined as follows:
$$
\begin{equation}
\psi_{e}=g(d)\psi_{e}^{+}+\psi_{e}^{-}
\label{eq:psi-e}
\tag{3}
\end{equation}
$$
where
$$
\begin{equation}
\psi_{e}^{+}=\frac{\lambda}{2}\langle\mathrm{tr}(\mathbf{\varepsilon})\rangle\_{+}^{2}+\mu \mathrm{tr}[\mathbf{\varepsilon}_{+}^{2}]
\label{eq:psi-e-positive}
\tag{4}
\end{equation}
$$

and
$$
\begin{equation}
\psi_{e}^{-}=\frac{\lambda}{2}\langle\mathrm{tr}(\mathbf{\varepsilon})\rangle\_{-}^{2}+\mu \mathrm{tr}[\mathbf{\varepsilon}_{-}^{2}]
\label{eq:psi-e-negative}
\tag{5}
\end{equation}
$$

Here, $\lambda$ and $\mu$ are the lame constant and shear moduli, respectively. The positive and negative bracket operators are given as follows:
$$
\begin{equation}
\langle x\rangle_{+}=\frac{x+|x|}{2}\\ ,\quad
\langle x\rangle_{-}=\frac{x-|x|}{2}
\label{eq:bracket}
\tag{6}
\end{equation}
$$

Thus the stress can be defined as follows:
$$
\begin{equation}
\mathbf{\sigma}=\frac{\partial\psi}{\partial\mathbf{\varepsilon}}
=g(d)[\lambda\langle\mathrm{tr}[\mathbf{\varepsilon}]\rangle_{+}\mathbf{I}+2\mu\mathbf{\varepsilon}\_{+}]+[\lambda\langle\mathrm{tr}[\mathbf{\varepsilon}]\rangle_{-}\mathbf{I}+2\mu\mathbf{\varepsilon}_{-}]
\label{eq:stress}
\tag{7}
\end{equation}
$$

The governing equaitons for this model are list below:
$$
\begin{equation}
\mathbf{\nabla}\cdot\mathbf{\sigma}=\mathbf{0}
\label{eq:stress-equilibrium}
\tag{8}
\end{equation}
$$
and
$$
\begin{equation}
\eta\frac{\partial d}{\partial t}=2(1-d)\mathcal{H}-\frac{\mathcal{G}_{c}}{l}(d-l^{2}\Delta d)
\label{eq:damage-equation}
\tag{9}
\end{equation}
$$
where $\eta$ is the viscosity coefficient. The history variable $\mathcal{H}$ is calculated as follows:

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
  "type":"msh4",
  "file":"tensile.msh",
  "savemesh":true
},
```
Here, you should use [gmsh](https://gmsh.info/) to generate the related msh file!

## Dofs
Next, you need to define the dofs, it should be $d$, $u_{x}$, and $u_{y}$ for 2d case, and $d$, $u_{x}$, $u_{y}$, and $u_{z}$ for 3d case. The `[dofs]` block should looks like:
```
"dofs":{
  "names":["d","ux","uy"]
},
```
or
```
"dofs":{
  "names":["d","ux","uy","uz"]
},
```
Here, the sequence of your `dofs` name matters !!!

## [elmts] and [mates]
The governings in Eq.$\eqref{eq:stress-equilibrium}$ and Eq.$\eqref{eq:damage-equation}$ can be implemented by using the following `elmts`:
```
  "elements":{
    "elmt1":{
      "type":"miehefracture",
      "dofs":["d","ux","uy"],
      "material":{
        "type":"miehefracture",
        "parameters":{
          "viscosity":1.0e-6,
          "Gc":2.7e-3,
          "eps":0.012,
          "K":121.15,
          "G":80.77,
          "stabilizer":1.0e-6,
          "finite-strain":false,
          "plane-strain":true
        }
      }
    }
  },
```
where `type=miehefracture` tells AsFem, our users want to call the phase-field fracture model.

For the material calculation, as mentioned before, several parameters are required, i.e., $E$, $\nu$ for the Youngs modulus and poisson ration, $\eta$, $\mathcal{G}_{c}$, and $l$ for the damage evolution.

<!-- It should be mentioned that, the de-coupled or staggered solution can be done easily by introducing the `usehist` parameter as follows:
```
[mates]
  [myfracmate]
    type=miehefracmate
    params=121.15 80.77 2.7e-3 0.012  1.0e-6     1
    //     lambda mu    Gc     L      viscosity  usehist
  [end]
[end]
```
then, your $\mathcal{H}$ will always use the $\mathcal{H}_{old}$ value from previous step. -->

## boundary condition
For the tensile test, one can use the following `bcs`:
```
"bcs":{
  "fixux":{
    "type":"dirichlet",
    "dofs":["ux"],
    "bcvalue":0.0,
    "side":["left","right"]
  },
  "fixuy":{
    "type":"dirichlet",
    "dofs":["uy"],
    "bcvalue":0.0,
    "side":["bottom"]
  },
  "loading":{
    "type":"dirichlet",
    "dofs":["uy"],
    "bcvalue":"0.1*t",
    "side":["top"]
  }
},
```
and for the shear failure test, one can use:
```
"bcs":{
  "fixux":{
    "type":"dirichlet",
    "dofs":["ux"],
    "bcvalue":0.0,
    "side":["bottom"]
  },
  "fixuy":{
    "type":"dirichlet",
    "dofs":["uy"],
    "bcvalue":0.0,
    "side":["left","right","bottom"]
  },
  "loading":{
    "type":"dirichlet",
    "dofs":["ux"],
    "bcvalue":"0.1*t",
    "side":["top"]
  }
},
```

Done? Yep, all the things is done!

# Run it in AsFem
Now, let's try your first phase-field fracture model in AsFem. You can create a new text file and name it as tensile.json or whatever you like. Then copy the following lines into your input file:
```
{
  "mesh":{
    "type":"msh4",
    "file":"tensile.msh",
    "savemesh":true
  },
  "dofs":{
    "names":["d","ux","uy"]
  },
  "elements":{
    "elmt1":{
      "type":"miehefracture",
      "dofs":["d","ux","uy"],
      "material":{
        "type":"miehefracture",
        "parameters":{
          "viscosity":1.0e-6,
          "Gc":2.7e-3,
          "eps":0.012,
          "K":121.15,
          "G":80.77,
          "stabilizer":1.0e-6,
          "finite-strain":false,
          "plane-strain":true
        }
      }
    }
  },
  "projection":{
    "type":"leastsquare",
    "scalarmate":["vonMises-stress"],
    "rank2mate":["stress","strain"]
  },
  "bcs":{
    "fixux":{
      "type":"dirichlet",
      "dofs":["ux"],
      "bcvalue":0.0,
      "side":["left","right"]
    },
    "fixuy":{
      "type":"dirichlet",
      "dofs":["uy"],
      "bcvalue":0.0,
      "side":["bottom"]
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
    "restarts":1600,
    "tolerance":1.0e-28
  },
  "nlsolver":{
    "type":"asfem",
    "maxiters":15,
    "abs-tolerance":5.0e-7,
    "rel-tolerance":1.0e-13,
    "s-tolerance":0.0
  },
  "timestepping":{
    "type":"be",
    "dt0":1.0e-4,
    "dtmax":5.0e-3,
    "dtmin":1.0e-12,
    "optimize-iters":5,
    "end-time":1.5e0,
    "growth-factor":1.05,
    "cutback-factor":0.75,
    "adaptive":true
  },
  "output":{
    "type":"vtu",
    "interval":5
  },
  "qpoints":{
    "bulk":{
      "type":"gauss-legendre",
      "order":2
    }
  },
  "postprocess":{
    "ux":{
      "type":"sideaveragevalue",
      "dof":"ux",
      "side":["top"]
    },
    "uy":{
      "type":"sideaveragevalue",
      "dof":"uy",
      "side":["top"]
    },
    "fx":{
      "type":"sideintegralrank2mate",
      "side":["top"],
      "parameters":{
        "rank2mate":"stress",
        "i-index":1,
        "j-index":1
      }
    },
    "fxy":{
      "type":"sideintegralrank2mate",
      "side":["top"],
      "parameters":{
        "rank2mate":"stress",
        "i-index":1,
        "j-index":2
      }
    },
    "fy":{
      "type":"sideintegralrank2mate",
      "side":["top"],
      "parameters":{
        "rank2mate":"stress",
        "i-index":2,
        "j-index":2
      }
    },
    "damage":{
      "type":"volumeaveragevalue",
      "dof":"d",
      "domain":["alldomain"]
    }
  },
  "job":{
    "type":"transient",
    "print":"dep",
    "restart":true
  }
}
```
You can also find the complete input file in `examples/pffracture/MieheTensile2d.json`.



If everything goes well, you can see the following image in your [Paraview](https://www.paraview.org/download/):
![shear](shear.jpeg)

![tensile](tensile.jpeg)
