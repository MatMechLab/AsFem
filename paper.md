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
date: 26 June 2022
bibliography: paper.bib

---

# Summary

The kinetics of multiphysics are naturally involved in the application of energy materials, such as lithium-ion batteries (LIBs) for example. Therefore, multiphysics modeling plays a crucial role in understanding the interactions between different fields in these applications, such as the chemo-mechanical interaction [@bai2021chemo;@bai2020chemo] and stress-enhanced species diffusion and its impact on cell performance[@bai2019two] in LIBs. The field of "multiphysics" or the "chemo-mechanics", aims to model the multiphysical problem in material sciences, for example, the diffusion of species, the phase transformation, and the large deformation under different fields. The major challenge of these models requires efficient numerical tools, easy-to-use parallelization, as well as a simple way to define the physical model. For this purpose, `AsFem`, a simple to use simulation package based on the finite element method [@ciarlet2002finite,@reddy2019introduction] for phase-field modeling and multiphysics coupling has been developed.


# Statement of need

`AsFem` is an open-source c++ package powered by the high-performance scientific computing library PETSc [@petsc-web-page]. By using C++, `AsFem` implement the object-oriented design for the finite element calculation without losing the flexibility and ease-of-use in parallelization. Due to the complexity of the finite element simulation procedure, it is important to simplify the model definition as well as the pre/post-processes. For this purpose, `AsFem` is intended to offer users the easy-to-use mathematic operators, i.e. the Laplacian operator and the inner product for vectors, as well as the tensor calculation (both the second-order tensor `RankTwoTensor class` and fourth-order tensor `RankFourTensor class`) operators and functions to simplify the model definition. Furthermore, instead of re-writing the code, users can easily switch their model from 2D to 3D by simply changing the mesh in the input file.

Why a new package is required? Indeed there are many open-source packages available, including FEniCS [@alnaes2015fenics], MFEM [@anderson2021mfem], deal.II [@arndt2021deal], and others. However, in those packages, users have to either offer weak forms (in FEniCS), or even the implementation of the elemental loop and the integration point loop (in deal.II). More importantly, users have to define the model by themselves. 

Therefore, with the emphasis on 'simple', AsFem has implemented the commonly used standard models for the chemo-mechanical phase-field modeling. Currently, `AsFem` has developed the finite strain model for hyperelastic materials (`MechanicsElmt`), the phase-field fracture model for brittle materials [@miehe2010phase] (`MieheFractureElmt`), the Cahn-Hilliard equation [@cahn1958free] and mechanically coupled spinodal decomposition [@bai2021chemo] (`CahnHilliardElmt` and `MechanicsCahnHilliardElmt`), the dendrite growth model [@kobayashi1993modeling] (`KobayashiElmt`), J2 plasticity model (`MechanicsElmt` with `J2PlasticityMaterial`), as well as the chemo-mechanical phase-field fracture model (`AllenCahnFractureElmt`), and more. Therefore, one can easily set up the simulation for the above-mentioned problems by simply modifying the parameters in the related input file and the mesh.

To further simplify the code development, AsFem offers the interface to different materials (`Materials`) for each `Element` (the PDEs/ODEs of your model). For instance, one can use the `MechanicsElmt` for all the solid mechanics problems together with different `Materials` for different analyses. In such a way, one can easily define the model by using specific PDEs/ODEs in AsFem with the different `Materials` for the calculation of different constitutive laws.

Moreover, users don't need to care about the choice of different basis functions or mesh for their model. Currently, AsFem offers the basis function for the interpolation up to second order, which can automatically detect the mesh type based on the mesh file. It supports the `Triangle` mesh with 3 (1st order) and 6 (2nd order) nodes, and the `Quadrilateral` mesh with 4 (1st order), 8, and 9 (2nd order) nodes in the 2D case. For the 3D case, the `Tetrahedron` mesh with 4 (1st order) and 10 (2nd order) nodes, and `Hexahedron` mesh with 8 (1st order) nodes, and 20 and 27 (2nd order) nodes are supported. Therefore, one can easily carry out the simulation on different mesh without the need for code modification. In most cases, users only need to modify the parameters in the input file for different models.

Furthermore, AsFem is not limited to chemo-mechanical phase-field equations, it implements the general user-defined-element (UEL) system (`User1Elmt` to `User20Elmt`) for the model definition, where one can easily write out the weak form and the related jacobian calculation on a single integration point (no need to write element/node loop and vector/matrix assemble). Accordingly, the user-defined-material (UMAT) system (`User1Material` to `User10Material`) has also been developed. Therefore, users can define their own model by using UEL and UMAT systems.

In such a way, `AsFem` enables the easy and rapid development of the physical model for FEM simulation. 


Part of `AsFem` features are list below:

- user-defined boundary condition (UBC)
- user-defined initial condition (UIC)
- user-defined element (UEL)
- user-defined material (UMAT)
- built-in 1st order and 2nd mesh generation for 1D, 2D, and 3D regular domain
- external mesh import, i.e., gmsh, netgen
- built-in postprocess
- vtu file output for results
- etc



# References 
