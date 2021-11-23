[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4603558.svg?style=flat-square)](https://doi.org/10.5281/zenodo.4603557)

# AsFem
A **simple** finite element method program, which is short for **AsFem**. AsFem is written in C++ and designed for phase-field modeling and multiphysics coupling. The [PETSc](https://www.mcs.anl.gov/petsc/) library is involved in AsFem for the efficient computing.

For '**simple**', we try to make the finite element programming and modeling, as simple as possible.

# Download
```
git clone https://github.com/M3Group/AsFem.git
```
For the detailed usage, one is referred to [AsFem Page](https://yangbai90.github.io/AsFem/) .

If one has access issues, the alternative link could be [AsFem-Gitee](https://gitee.com/m3group/AsFem.git), the git clone can be done via:
```
git clone https://gitee.com/m3group/AsFem.git
```


## Installation
The installation details of AsFem can be found here [AsFem-Installation](https://yangbai90.github.io/AsFem/install) .

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

# Tutorial
The tutorial is avialable here https://yangbai90.github.io/AsFem/Tutorial/step-0/

For Chinese users, the video lecture is available on bilibili, please see [AsFem-Lecture](https://space.bilibili.com/100272198/channel/detail?cid=193605).

# Document
The code is documented by the [Doxygen](https://www.doxygen.nl/index.html) package, one can generate the pdf file or html files via:
```
doxygen
```

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

# Discussion & bug report
If you are interested in AsFem or have any questions, just feel free to send me an email [Mail2Me](mailto:yangbai90@outlook.com) or join the QQ group for more discussion .
```
QQ group: 879908352
```
