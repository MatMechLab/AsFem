---
title: AsFem Release
date: 2020-04-19 16:20:52
---

All the release can be download from here [AsFem-Release](https://github.com/yangbai90/AsFem/releases)

## 💥 Breaking Changes

- Shift the whole framework to [PETSc](https://www.mcs.anl.gov/petsc/) to achieve the better nonlinear solver and more efficient parallelization
- Add support for results ouput with folder choice
- Remove the linear solver from [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## ⭐ Features

- Linear/Nonlinear elastic problem analysis
- Phase-field fracture model in 2D/3D
- Dendrite growth model in 2D
- Thermo-mechanical coupling in 2D/3D

## 🐞 Bug Fixes

- Fix the incorrect order of gauss integration when using triangle-3 node from Gmsh

## 🛠 Improvements
- Better nonlinear equation solvability from SNES in PETSc
- Improved parallelization with mumps solver
