---
title: AsFem manual
date: 2020-03-31 20:22:15
tags:
- document
categories:
- document

mathjax: true
---
# Before we start
Currently, some available demos can be found here:

https://space.bilibili.com/100272198/channel/detail?cid=90241
 
## Write a simple input file:

structure of the block in your input file should look like:
```
[block_name]
type=...
option1=...
option2=...
[end]
```
Here we use the [blockname]/[end] bracket pair to set the properties we want. Thus each block pair must end up with an '[end]'. Remember, even there are sub-blocks inside one big block, each block itself still must end up with an '[end]'. Otherwise, your input file will complain about errors to you!

For one minimal input file, you need:
```
[mesh]
...
[end]
```

and the name or name list of the defined degrees of freedom (DoFs), which will be displayed in the Paraview when you open the result file(vtu file). Besides, the name of related DoFs is required when you want to apply the boundary condition.
```
[dofs]
...
[end]
```

and the element or the module you want to use:
```
[elmts]
...
[end]
```

and the information for the analysis(we call it the [job] block)
```
[job]
...
[end]
```

In summary, you need at least:
```
[mesh]
...
[end]
[dofs]
...
[end]
[elmts]
...
[end]
[job]
...
[end]
```

<!--###########################################################################-->
<!--For the first chapter-->
<!--###########################################################################-->

# Installation

## Windows
For windows users, it is recommended to download the executable file from here [AsFem-Win](https://github.com/yangbai90/AsFem/releases)

If someone really needs to compile the source code by themselves, one can do the following steps:

### step-1
Firstly, you need to install the PETSc in your windows system. In order to do this, you need to install [Cygwin](http://www.cygwin.com/). Download it and double click on it, it will ask you to select the FTP server(please choose the one close to you to speed up the download rate).

After that, you can search the GCC, g++, python, make components (GCC/g++, make, python are required to install PETSc). After the necessary components are selected, Cygwin will download it and install them into the folder you selected. If everything is fine, you should find the 'mintty.exe' file in your Cygwin-install-path/bin folder.

### step-2
In order to let the Cygwin compiler knows that we have the msvc c/c++ compiler, please open your start->visualstudo2019->x86_64-native-compiler (whatever the name it is, just a terminal command-prompt, open it!). After the terminal is open, cd into the 'mintty.exe' file's folder and run
```
/usr/bin/bash --login
```
it will open a new terminal (named by Cygwin maybe), now both the win32fe compiler and GCC compiler are ready for the next step. Remember, in the following steps, you should do all the compilation stuff in this new terminal!!!

### step-3
Now we can do the installation for PETSc in windows since your new terminal is a Linux-type terminal, you can download the PETSc via:
```
curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.4.tar.gz
```
or download PETSc from your favorite browser, and choose the version you like from here: https://www.mcs.anl.gov/petsc/mirror/release-snapshots/

Then you can unzip it via 'tar -xf petsc.tar.gz' or just right-click on it and unzip it. After that, cd into the petsc folder from your new terminal, and run the configuration file like below:
```
./configure --with-cc='win32fe cl' --with-fc=0 --with-cxx='win32fe cl' --with-mpi=0 \
--download-f2cblaslapack \
CFLAGS='-fPIC -fopenmp -O3 -march=native' \
CXXFLAGS='-fPIC -fopenmp -O3 -march=native' \
FFLAGS='-fPIC -fopenmp -O3 -march=native' \
FCFLAGS='-fPIC -fopenmp ' \
F90FLAGS='-fPIC -fopenmp ' \
F77FLAGS='-fPIC -fopenmp '
```
after the configuration is done, it will ask you to do several 'make xxxxxx' steps, then you just need to copy the 'make xxx' command from the terminal and copy it and run it.
After several steps (3~4, normally it is 1) 'make xxxx' 2) 'make xxx install' 3) 'make xxxx test'), your PETSc should be ready. By default, the petsc is installed in the same folder as your petsc source code, but it may have a different name, i.e. win-opt-arch, something like that.

### step-4
Now you can compile AsFem after you download the source code, you need to edit the CMakeLists.txt file in AsFem folder, please choose the windows release configuration in CMakeLists.txt file, and make sure the path to your PETSc and Eigen is correct (in windows, the path is 'C:\\Program\\xxx', however, it should be modified as 'C:/Program/XXX', otherwise it can not be recognized). Another important thing is, in Linux the library file is named as 'xxx.so'. However, it is different in windows, where you should give 'petsc.lib' or whatever the name you see in your petsc's lib folder.

After your CMakeLists.txt file is ready, open your VisualStudio2019 and open the AsFem folder, then compile it.
In order to run asfem, you need to put the petsc.dll file together with asfem. 

## Linux

For Linux users, it is much simpler. All you need is the g++ compiler and one working mpi compiler. Then you can do:

### step-1
install the PETSc

### step-2
download the source code:
```
git clone https://github.com/yangbai90/AsFem.git
```
then go to the source code folder and edit the CMakeLists.txt file, change the path of PETSc and MPI to your own one. Then you can run:
```
cmake CMakeLists.txt && make -j8
```
after the compilation is finished, the executable file 'asfem' should be found under the bin folder.

<!--###########################################################################-->
<!--For the second chapter-->
<!--###########################################################################-->
# Mesh generation

For the mesh generation, the [mesh] block is required.
For example, if one wants to generate the 1d Lagrange mesh, you can use:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=10
  meshtype=edge2
[end]
```
if one wants to save the generated mesh, then the option 'savemesh=true' is required:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=10
  meshtype=edge2
  savemesh=true
[end]
```
The mesh will be saved as a *.vtu* file, which should have the name: your_input_file_name+'_mesh.vtu'(*.i* is removed from your input file name). For example, if your input file is: *test.i*, then the mesh file name is: *test_mesh.vtu*.

Additionally, AsFem will also print this message in your console:
```
**************************************************************
*** Start to read input file ...                           ***
***   start to crate mesh ...                              ***
***     save mesh to                   hex20m1_mesh.vtu    ***
***   mesh generation finished !                           ***
```

Moreover, if one wants to print out the mesh information in the terminal, one can use the 'printmesh=true' option. For example, if we want to print the information of a 1D edge3 type Lagrange mesh, we can use:
```
[mesh]
  type=asfem
  dim=1
  xmin=0.0
  xmax=1.0
  nx=5
  meshtype=edge4
  savemesh=false
  printmesh=true
[end]
```
then the output should look like:
```
***-------------------------------------------------------------------***
*** Mesh information summary:                                         ***
***  Nodes=       16, Elmts=        7, NodesPerBulkElmt=  4           ***
***  Max dim= 1, Min dim= 0, PhyGroup=    3, Meshtype= edge4, Order=3 ***
***  Physical id                    Phsical Name             Elmts    ***
***        1                                left                 1    ***
***        2                               right                 1    ***
***        3                           alldomain                 5    ***
***-------------------------------------------------------------------***
```
if one wants to print out the details of the mesh, then 'printmesh=dep' option should be used. Consequently, the output should look like:
```
***-------------------------------------------------------------------***
*** Mesh information summary:                                         ***
***  Nodes=       16, Elmts=        7, NodesPerBulkElmt=  4           ***
***  Max dim= 1, Min dim= 0, PhyGroup=    3, Meshtype= edge4, Order=3 ***
***  Physical id                    Phsical Name             Elmts    ***
***        1                                left                 1    ***
***        2                               right                 1    ***
***        3                           alldomain                 5    ***
***  Physical group information (ID and element ID)                   ***
***  phyname=                     left, element id=        1          ***
***  phyname=                    right, element id=        2          ***
***  phyname=                alldomain, element id=        3          ***
***  phyname=                alldomain, element id=        4          ***
***  phyname=                alldomain, element id=        5          ***
***  phyname=                alldomain, element id=        6          ***
***  phyname=                alldomain, element id=        7          ***
***  Element connectivity information(element id: node index):        ***
***  elmt id=        1:       1                                       ***
***  elmt id=        2:      16                                       ***
***  elmt id=        3:       1        2        3        4            ***
***  elmt id=        4:       4        5        6        7            ***
***  elmt id=        5:       7        8        9       10            ***
***  elmt id=        6:      10       11       12       13            ***
***  elmt id=        7:      13       14       15       16            ***
***  Node coornidates (node id, x, y, z and weight)                   ***
***          1:   0.0000e+00,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          2:   6.6667e-02,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          3:   1.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          4:   2.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          5:   2.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          6:   3.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          7:   4.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          8:   4.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***          9:   5.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         10:   6.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         11:   6.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         12:   7.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         13:   8.0000e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         14:   8.6667e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         15:   9.3333e-01,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***         16:   1.0000e+00,   0.0000e+00,   0.0000e+00,   1.000e+00 ***
***-------------------------------------------------------------------***
```
Since **printmesh=dep** option will print out all the information, includes the nodal coordinates and elemental connectivity, it is suggested **not** use this option in your [mesh] block. If one wants to check the mesh, 'savemesh=true' is already good enough!

Similarly, for 2D and 3D cases, one can use:
```
[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=10
  ny=10
  meshtype=quad4
[end]
```
and
```
[mesh]
  type=asfem
  dim=3
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=1.0
  nx=10
  ny=10
  nz=10
  meshtype=hex8
[end]
```
Or, one can also use:
```
[mesh]
  type=asfem
  dim=3
  nx=10
  ny=10
  nz=10
  meshtype=hex8
[end]
```
then the unit([0,1] in 1D, [0,1]x[0,1] in 2D, [0,1]x[0,1]x[0,1] in 3D) domain will be used by default.

Currently, AsFem supports the following kinds of built-in mesh:

- 1D--> edge2, edge3, edge4

- 2D--> quad4, quad8, quad9

- 3D--> hex8, hex20, hex27

<!--###########################################################################-->
<!--For Degrees of freedom (DoFs)-->
<!--###########################################################################-->
# Degrees of freedom (DoFs)
In AsFem, the name of DoFs is required before the simulation can start. Once the name of DoFs is given, it is extremely easy to apply the boundary condition and also define the elements we want to use.
The [dofs] block is very simple, the complete block should look like:
```
[dofs]
name=ux uy uz
[end]
```
this block tells AsFem that, we defined 3 DoFs, their names are: ux(the first one), uy(the second one), uz(the third one). Thus each node of our mesh will have 3 DoFs. 

By the way, when the result file is written to a VTU file, the related result will be associated with the name you defined. Thus it is very convenient to do the postprocessing.

<!--###########################################################################-->
<!--For Elements-->
<!--###########################################################################-->
# Elements(uel)

The available elements in AsFem are implemented as the user-defined element (uel) way. Currently, it supports
- [X] Poisson Equation (type=poisson)
- [x] Solid Mechanics (2D and 3D) (type=mechanics)
- [x] Fick Diffusion (type=diffusion)
- [x] Wave Propagation (type=wave)
- [x] Solid Mechanics (2D and 3D) (type=mechanics)
- [x] Cahn Hilliard diffusion (type=cahnhilliard)
- [x] Mechanically Coupled Thermal Analysis (type=thermalmechanics)
- [ ] Phase Field Fracture (under developing, type=linearelasticphasefieldfracture)

The related block in your input file should look like:
```
[elmts]
  [your_sub_block_name]
    type=your_element_name
    dofs=ux uy
    mate=name_of_your_material  // this can be ignored
    domain=alldomain            // domain-id, can be ignored
  [end]
[end]
```
For example, if we want to the linear elastic analysis, one can use:
```
[elmts]
  [solids]
    type=mechanics
    dofs=ux uy
    mate=elastic
    domain=alldomain
  [end]
[end]
```
then the related material block can be defined as following:
```
[mates]
  [elastic]
    type=linearelastic
    params=1.0e5 0.25
    // params= E nu(Youngs modulus and poisson ratio)
  [end]
[end]
```
where, 'type=' should be given firstly, and the type name should be the buil-in type name or user defined one, i.e. user1, user2... 

Once the type is given, the 'dofs=' should be list, otherwise, nobody knows what you want to do. Besides, if one wants to define the material, then 'mate=mate_block_name' should be given. 

Remember, the 'mate=' is given by **the name of the material block name**, _not the material type name in [mates] block_.

In [elmts] block, the domain name and material name can be ignored, however, if one wants to use a different material or your own material model, then a [mates] block is required to define the material.

<!--###########################################################################-->
<!--For Materials(umat)-->
<!--###########################################################################-->
# Materials(umat)

It is general and also convenient to has an interface to allow users to define their own material model, i.e. the constitutive law for solid mechanics or the free energies for phase-field simulation, hereby, the [mates] block is introduced in AsFem.

The structure of the [mates] block looks like:
```
[mates]
  [mate1]
    type=mate_name1
    params=val1 val2 val3
  [end]
  [mate2]
    type=mate_name2
    params=x1 y1 z1
  [end]
[end]
```
where the user can define multiple sub-blocks ([mate1],[mate2])inside the [mates]/[end] block pair, the **sub block name should be different**, otherwise AsFem will complain about the duplication. The 'type=' indicates the material you want to use, while 'params=' allows the user to define the parameters they want to use for this material.

Currently, AsFem supports:
- [ ] Constant Poisson Material (type=constpoisson)
- [ ] Nonlinear Poisson Material (type=nlpoisson)
- [ ] Tensor Poisson Material (type=tensorpoisson)
- [ ] Constant Diffusion Material (type=constdiffusion)
- [ ] Nonlinear Diffusion Material (type=nldiffusion)
- [ ] Tensor Diffusion Material (type=tensordiffusion)
- [ ] Constant Wave Propagation Material (type=constwave)
- [ ] Linear Elastic Material (type=linearelastic)
- [ ] Saint Venant Elastic Material (type=saintvenant)
- [ ] Neo-Hookean Elastic Material (type=neohookean)
- [ ] Mooney-Rivlin Elastic Material (type=mooneyrivlin)
- [ ] Wave Propagation Material (type=constwave)
- [ ] Thermal Mechanics Coupling Material (type=thermelastic)
- [ ] Linear Elastic Phase Field Fracture Material (type=lpffracture)
- [ ] Linear Elastic Cahn-Hilliard Diffusion Material (type=elasticch)
- [ ] User Defined Material (UMAT,currently AsFem suppor 20 umat) (type=user1[,user2,user3,...,user20])

In order to reduce the difficulty of coding, AsFem introduces the tensor calculation, which is quite convenient for tensor-based constitutive law and other tensor calculation.


<!--###########################################################################-->
<!--For Boundary conditions(BC)-->
<!--###########################################################################-->
# Boundary conditions(BCs)
Boundary conditions are quite important if someone has to have a numerical solution. For most of the PDEs, the boundary condition is the key point to let you have the correct solution. In this chapter, the boundary condition system of AsFem will be introduced.

Normally, there are three types of boundary conditions: 1) Dirichlet boundary conditions, where the degree of freedom (Dof) equals to some specific values at some specific boundary, namely $$u_{i}=u_{0}$$ . 2) Neumann boundary conditions, where the gradient of the degree of freedom equals to some certain values at some specific boundary, namely $\nabla u\cdot\vec{n}=g_{0}$ . 3) Robin boundary conditions, where the combination of the previous two types is introduced, namely $\nabla u\cdot\vec{n}=f(u)$ , where $f(u)$ is the function of our DoF, it could be linear or nonlinear.

Thus, one can apply the different boundary conditions in AsFem as follow:
```
[bcs]
  [mybc1]
    type=dirichlet[neumann]
    dof=ux
    boundary=bottom
    value=0.0
  [end]
[end]
```
where the [bcs]/[end] pair denotes the boundary condition block, inside this master block, one can define several sub-blocks for different boundary conditions. The basic information for a sub-block looks like:
```
[mybc1]
  type=dirichlet[neumann,user1,user2,...]
  dof=ux
  boundary=bottom
  value=0.0
[end]
```
where 'type=' defines which type of boundary condition one wants to use, currently dirichlet, neumann, user1~user10 is supported. 

For user1-user10, they are designed to allow users to do any kind of boundary conditions they want, namely the user-defined boundary condition, short for ubc. Thus, if one wants to apply a robin boundary condition, then one can use one of the user1-user10 cpp file in src/BCs folder to implement your own BC.

The 'dof=' defines the name of the dof for the applied boundary condition. **The dof name should be one of the name in [dofs]/[end] block**. 

While 'boundary=' represents the name of the boundaries we want to apply, yes, the **boundaries** not the boundary. It means, if you want to set 'u=1.0' at both the left, right, bottom and top edges of a rectangle domain, then one can use 'boundary=left right bottom top'. 

The last one, 'value=' is the boundary value one wants to apply.

## Dynamic loading
So, what should I do if I want to apply a dynamic displacement loading to my cubic box? The answer is, **'t'**. Yes, just put the time-dependent expression, i.e. 0.5*t in your 'value=', AsFem will treat it as a time-dependent boundary condition. For sure, for this case, you definitely need a transient analysis, not a static one!

For example, the following lines defined a dynamic displacement loading (ux) on the top edge of a rectangle domain, but fixed other edges:
```
[bcs]
  [fixUx]
    type=dirichlet
    dof=ux
    boundary=bottom
    value=0.0
  [end]
  [fixUy]
    type=dirichlet
    dof=uy
    boundary=bottom left right top
    value=0.0
  [end]
  [loadUx]
    type=dirichlet
    dof=ux
    value=1.0*t
    boundary=top
  [end]
[end]
```


<!--###########################################################################-->
<!--For Initial conditions (ICs)-->
<!--###########################################################################-->
# Initial conditions (ICs)

For time-dependent problems, the initial condition (IC) is very important to have a good result. Similar to the boundary conditions, the initial conditions in AsFem is applied by the [ics]/[end] block, which looks like:
```
[ics]
  [myic1]
    type=const
    dof=d
    params=0.0
  [end]
[end]
```
where the sub-block inside the [ics]/[end] master pair is used to define multiple initial conditions. The basic information of a initial condition is defined as:
```
[myic1]
  type=const
  dof=d
  params=0.0
  domain=alldomain
[end]
```
where 'type=' denotes which type of initial condition you want to use. Currently, AsFem supports 'const, random, circle, rectangle, user1~user10' initial conditions.

The 'dof=dof_name' defines the dof we want to apply the IC, while 'params=' is the parameter list for this IC.

The 'domain=' denotes, which domain we want the current IC to be applied. **Once 'domain=' is ignored, all the domain will be used !!!**

For example, if we want to let c has a randomly distribution in [0.6,0.66], while mu to be the constant (0.0), we can use:
```
[ics]
  [myc]
    type=random
    dof=c
    params=0.6 0.66
  [end]
  [mymu]
    type=const
    dof=mu
    params=0.0
  [end]
[end]
```

<!--###########################################################################-->
<!--For Examples -->
<!--###########################################################################-->


That's all, enjoy it!