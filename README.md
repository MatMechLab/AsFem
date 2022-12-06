[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4603558.svg?style=flat-square)](https://doi.org/10.5281/zenodo.4603557)

# AsFem

A **simple** finite element method program, which is short for **AsFem**. AsFem is written in C++ and designed for phase-field modeling and multiphysics coupling. The [PETSc](https://www.mcs.anl.gov/petsc/) library is involved in AsFem for the efficient computing.


# Download

```
git clone https://github.com/M3Group/AsFem.git
```
For the detailed usage, one is referred to [AsFem Page](https://m3group.github.io/AsFem/) .

If one has access issues, the alternative link could be [AsFem-Gitee](https://gitee.com/m3group/AsFem.git), the git clone can be done via:
```
git clone https://gitee.com/m3group/AsFem.git
```


# Installation

The installation details of AsFem can be found here [AsFem-Installation](https://m3group.github.io/AsFem/install/) .

After you've installed your PETSc (PETSC_DIR) and MPI (MPI_DIR), all you have to do is:
```
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
make -j4
```

# Demos
AsFem has several demos, you can visit them either on [Bilibili](https://space.bilibili.com/100272198/channel/detail?cid=90241) or on [Youtube](https://www.youtube.com/playlist?list=PLVEpIo_wvYmaLPoLjj5Lg93YvYy9flkN8).

spinodal-decomposition              |  double-notch failure
:-------------------------:|:-------------------------:
![](figures/CahnHilliard.gif)      |  ![](figures/DoubleNotch.gif)

# Tutorial & examples
The tutorial is avialable here https://m3group.github.io/AsFem/Tutorial/step-0/

For Chinese users, the video lecture is available on bilibili, please see [AsFem-Lecture](https://space.bilibili.com/100272198/channel/detail?cid=193605).

Currently, one can find the tested input files in both the `AsFem/examples` folder as well as `AsFem/test_input` folder. For a better understanding of how to use AsFem, please take a look at these input files first.

# Document
The code is documented by the [Doxygen](https://www.doxygen.nl/index.html) package, one can generate the pdf file or html files via:
```
doxygen
```

More documents and details can be found in AsFem's [homepage](https://m3group.github.io/AsFem/).

# Author
[Yang Bai](https://yangbai90.github.io/)


# Citation
```
@misc{yang_bai_2021_4603558,
  author       = {Yang Bai},
  title        = {{AsFem: a simple finite element method program for
                   phase-field modeling and multiphysics coupling}},
  month        = feb,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {master},
  doi          = {10.5281/zenodo.4603558},
  url          = {https://doi.org/10.5281/zenodo.4603558}
}
```

# Contribute

If you discover bugs in the `AsFem` package, please create an issue at the project repository on GitHub at https://github.com/M3Group/AsFem.

If you find the AsFem package useful, we welcome your code and documentation contributions. To contribute, fork the repository on GitHub, and submit a pull request at https://github.com/M3Group/AsFem.


# Contact

Submit bug reports and feature requests on the [repository issue tracker](https://github.com/M3Group/AsFem/issues).


# Discussion

If you are interested in AsFem or have any questions, just feel free to send me an email [Mail2Me](mailto:yangbai90@outlook.com) or join the QQ group for more discussion .
```
QQ group: 879908352
```
