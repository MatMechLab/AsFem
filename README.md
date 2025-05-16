[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4603558.svg)](https://doi.org/10.5281/zenodo.4603557)

# AsFem

**AsFem**-**A**dvanced **S**imulation kit based on **F**inite **E**lement **M**ethod. AsFem is written in C++ and designed for phase-field modeling and multiphysics coupling. The [PETSc](https://www.mcs.anl.gov/petsc/) library and [MPI](https://www.open-mpi.org/) package are involved in AsFem for the scalable computing.


# Download

```
git clone https://github.com/MatMechLab/AsFem
```
For comprehensive details on usage, please refer to the [AsFem Page](https://matmechlab.github.io/AsFem/).

If you encounter any access issues, an alternative link to explore is [AsFem-Gitee](https://gitee.com/matmechlab/AsFem). Git cloning can be executed using the following command:
```
git clone https://gitee.com/matmechlab/AsFem.git
```


# Installation

Detailed installation instructions for AsFem are available at [AsFem-Installation](https://matmechlab.github.io/AsFem/install/).

Once you have installed PETSc and MPI (set up your PETSC_DIR and MPI_DIR), the next steps are as follows:
```
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
make -j4
```

# Demos

You can explore several AsFem demos on platforms like [Bilibili](https://space.bilibili.com/100272198/channel/collectiondetail?sid=73185) and [Youtube](https://www.youtube.com/playlist?list=PLVEpIo_wvYmaLPoLjj5Lg93YvYy9flkN8).


spinodal-decomposition              |  double-notch failure
:-------------------------:|:-------------------------:
![](figures/CahnHilliard.gif)      |  ![](figures/DoubleNotch.gif)

# Tutorial & examples

The tutorial can be accessed at this link: [AsFem Tutorial](https://matmechlab.github.io/AsFem/Tutorial/step-0/)

Chinese users can access the video lecture on Bilibili at [AsFem-Lecture](https://space.bilibili.com/100272198/channel/collectiondetail?sid=73182).

At present, the tested input files for AsFem are located in two directories: `AsFem/examples` and `AsFem/test_input`. For an enhanced understanding of how to utilize AsFem effectively, it is recommended to first examine these input files.

# Document

The code is documented using the [Doxygen](https://www.doxygen.nl/index.html) package. You can generate PDF or HTML files using the following commands:
```
doxygen
```

Additional documentation and details about AsFem are available on its [homepage](https://matmechlab.github.io/AsFem/)

# Author
[Yang Bai](https://yangbai90.github.io/)

[MMLab members](https://www.x-mol.com/groups/matmechlab/people)

# Funding
This project is supported by:

- Start up fund of the Great Bay University
- General Program of Guangdong Natural Science Foundation, 2025-2027
- Youth Program of the National Natural Science Foundation of China, 2025-2027


# Citation
```
@misc{yang_bai_2021_4603558,
  author       = {Yang Bai},
  title        = {{AsFem: Advanced Simulation kit based on Finite 
                   Element Method}},
  month        = mar,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {V0.8},
  doi          = {10.5281/zenodo.4603558},
  url          = {https://doi.org/10.5281/zenodo.4603558}
}
```

# Contribute

If you discover bugs in the `AsFem` package, please create an issue at the project repository on GitHub at https://github.com/MatMechLab/AsFem.

If you find the AsFem package useful, we welcome your code and documentation contributions. To contribute, fork the repository on GitHub, and submit a pull request at https://github.com/MatMechLab/AsFem.


# Contact

Submit bug reports and feature requests on the [repository issue tracker](https://github.com/MatMechLab/AsFem/issues).


# Discussion

If you have an interest in AsFem or any questions, please don't hesitate to [Mail2Me](mailto:yangbai90@outlook.com), or join our QQ group for further discussions.
```
QQ group: 879908352
```
