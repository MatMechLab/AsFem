---
title: 'AsFem: A simple to use finite element package for phase-field modeling and multiphysics simulation'
tags:
  - c++
  - phase-field
  - finite element method
  - simulation
  - solid mechanics
  - computational materials
authors:
  - name: Yang Bai^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-1946-4377
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Yunxia Li #
    affiliation: "2"
affiliations:
 - name: Department of Microstructure Physics and Alloy Design, Max-Planck-Institut für Eisenforschung GmbH, Max-Planck-Strasse 1, 40237 Düsseldorf, Germany
   index: 1
 - name: Faculty of Arts and Humanities, University of Cologne, Dürener Straße 56–60, 50931 Köln, Germany
   index: 2
date: 16 February 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The kinetics of multiphysics are naturally involved in the application of energy materials, such as lithium-ion batteries (LIBs) for example. Multiphysics modeling plays a crucial role in understanding the interactions between fields in these applications, such as the chemo-mechanical interaction in LIBs [@bai2021chemo,@bai2020chemo]. The field of "multiphysics" or the "chemo-mechanics", aims to model the multiphysical problem in material science, for example, the diffusion of species, the phase transformation, and the large deformation under different fields. The major challenge of these models requires efficient numerical tools, easy-to-use parallelization, as well as a simple way to define the physical model. For this purpose, `AsFem`, a simple to use finite element method (FEM) simulation package for phase-field modeling and multiphysics coupling has been developed.


# Statement of need

`AsFem` is an open-source c++ package powered by the high-performance scientific computing library PETSc [@petsc-web-page]. By using C++, `AsFem` implement the object-oriented design for the finite element calculation without losing the flexibility and ease-of-use in parallelization. Due to the complexity of the finite element simulation procedure, it is important to simplify the model definition as well as the pre/post-processes. For this purpose, `AsFem` is intended to offer users the easy-to-use mathematic operators, i.e. the Laplacian operator and the inner product for vectors, as well as the tensor calculation (both the second-order tensor and fourth-order tensor) operators and functions to simplify the model definition. Furthermore, `AsFem` offers the user-defined-element (UEL) and user-defined-material (UMAT) to further minimize complexity, where all calculations are centered on the integration point (gauss point). Therefore, there is no need to write the nodal loop or the elemental loop for their calculation, which is still a standard procedure in other packages. Furthermore, instead of re-writting the code, users can easily swap their model from 2D to 3D by simply changing the mesh in the input file.

Currently, `AsFem` has implemented the phase-field fracture model for brittle materials [@miehe2010phase], mechanically coupled spinodal decomposition [@bai2021chemo], the dendrite growth model [@kobayashi1993modeling], J2 plasticity model, as well as the chemo-mechanical phase-field fracture model, and more. Besides that, users can easily implement their own model by either combing the built-in models or using the UEL/UMAT system. In such a way, `AsFem` will enable the easy and rapid development of the physical model for FEM simulation. 


Part of `AsFem` features are list below:

- user-defined boundary condition (ubc)
- user-defined initial condition (uic)
- user-defined element (UEL)
- user-defined material (UMAT)
- built-in 1st order and 2nd mesh generation for 1D, 2D, and 3D regular domain
- external mesh import, i.e., gmsh, netgen
- built-in postprocess
- vtu file output for results
- etc



# References 
