---
title: Installation
date: 2020-04-15 15:49:58
categories:
- Documents
tags:
- Documents
- Installation
---

# Windows
For windows users, it is highly recommended to download the executable file from here [AsFem-Win](https://github.com/yangbai90/AsFem/releases)

If someone really needs to compile the source code by themselves, one can do the following steps:

## step-1
Firstly, you need to install the PETSc in your windows system. In order to do this, you need to install [Cygwin](http://www.cygwin.com/). Download the exe file and double click on it, then you will be asked to select the FTP server(please choose the one close to you to speed up the download rate).

After that, you can search the GCC, g++, python, make components (GCC/g++, make, python are required to install PETSc). After the necessary components are selected, Cygwin will download it and install them into the folder you selected. If everything is fine, you should find the 'mintty.exe' file in your Cygwin-install-path/bin folder.

## step-2
In order to let the Cygwin compiler knows that we have the msvc c/c++ compiler(For sure, you must install the [VS2019](!https://visualstudio.microsoft.com/)), please open your start->visualstudo2019->x86_64-native-compiler (whatever the name it is, just a terminal command-prompt, open it!). Once the terminal is opened, cd into the 'mintty.exe' file's folder and run(please copy the following line into your terminal completely):
```
/usr/bin/bash --login
```
it will open a new terminal automatically (named by Cygwin maybe), now both the win32fe compiler and GCC compiler are ready for the next step. Remember, in the following steps, you should do all the compilation stuff in this new terminal!!!

## step-3
Now we can do the installation for PETSc in windows since your new terminal is already a Linux-type terminal, you can download the PETSc via:
```
curl -L -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.11.4.tar.gz
```
or download PETSc from your favorite browser, and choose the version you like from here: https://www.mcs.anl.gov/petsc/mirror/release-snapshots/

Then you can unzip it via ‘tar -xf petsc.tar.gz’ or just right-click on it and unzip it. After that, cd into the petsc folder from your new terminal, and run the configuration file like below(please don't download any external packages, otherwise the configuration will fail):
```
./configure --with-cc='win32fe cl' --with-fc=0 --with-cxx='win32fe cl' \
--with-mpi=0 \
--with-openmp=1 \
--with-debugging=0 \
--download-f2cblaslapack \
COPTFLAGS='-O3 -fopenmp' \
CXXOPTFLAGS='-O3 -fopenmp' \
PETSC_DIR=`pwd`
```
after the configuration is done, it will ask you to do several 'make xxxxxx' steps, then you just need to copy the 'make xxx' command from the terminal and copy it and run it.
After several steps (3~4, normally it is 1) 'make xxxx' 2) 'make xxx install' 3) 'make xxxx test'), your PETSc should be ready. By default, the petsc is installed in the same folder as your petsc source code, but it may have a different name, i.e. win-opt-arch, something like that.

## step-4
Now we can compile the AsFem source code. All you need to do is editing the CMakeLists.txt file in AsFem's source code folder. Please choose the windows release configuration in CMakeLists.txt file, and make sure the path to your PETSc and Eigen is correct (in windows, the path is 'C:\\Program\\xxx', however, it should be modified as 'C:/Program/XXX', otherwise it can not be recognized). Another important thing is, in Linux the library file is named as 'xxx.so'. However, it is different in windows, where you should give 'petsc.lib' or whatever the name you see in your petsc's lib folder.

After your CMakeLists.txt file is ready, open your VisualStudio2019 and open the AsFem folder(here you don't need to create a new project, just open the AsFem folder should be fine), then compile it.
In order to run asfem, you need to put the petsc.dll file together with asfem executable.


# Linux

For Linux users, it is much simpler. All you need is the g++ compiler and one working MPI compiler. Then you can do:

## step-1
install the PETSc
### step-1-1
run the configuration file for petsc:
```
./configure \
--prefix=path-to-your-petsc-installation-folder \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=1 \
--with-mpi=1 \
--with-shared-libraries=1 \
--with-cxx-dialect=C++14 \
--with-fortran-bindings=0 \
--with-sowing=0 \
--download-hypre=1 \
--download-fblaslapack=1 \
--download-metis=1 \
--download-ptscotch=1 \
--download-parmetis=1 \
--download-superlu_dist=1 \
--download-scalapack=1 \
--download-mumps=1 \
--download-slepc=1 \
COPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
CXXOPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
FOPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
PETSC_DIR=`pwd`
```
or you can also compile a minimal petsc(without any external packages):
```
./configure \
--prefix=path-to-your-petsc-installation-folder \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=1 \
--with-mpi=1 \
--with-shared-libraries=1 \
--with-cxx-dialect=C++14 \
--with-sowing=0 \
--download-fblaslapack=1 \
COPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
CXXOPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
FOPTFLAGS='-fPIC -fopenmp -O3 -march=native -mtune=native ' \
PETSC_DIR=`pwd`
```

if everthing is fine, you will be asked to run two to three times for compilation and installation, the command line should looks like:
for compilation:
```
make PETSC_DIR=/home/user/Program/temp/petsc-3.12.5 PETSC_ARCH=arch-linux-c-opt all
```
and for installation:
```
make PETSC_DIR=/home/user/Program/temp/petsc-3.12.5 PETSC_ARCH=arch-linux-c-opt install
```
then for the tests:
```
make PETSC_DIR=/home/user/Program/MooseLib/petsc/3.12.5 PETSC_ARCH="" test
```
Now the petsc is ready, we can compile the AsFem library.

## step-2
We can download the AsFem source code via:
```
git clone https://github.com/yangbai90/AsFem.git
```
then go to the source code folder and edit the CMakeLists.txt file, change the path of PETSc and MPI to your own one. Then you can run:
```
cmake CMakeLists.txt && make -j8
```
after the compilation is finished, the executable file 'asfem' should be found under the bin folder.

## step-3
Add the 'asfem' executable to your bashrc, which can allow you to run the asfem from terminal directly:
Open your bashrc:
```
nano ~/.bashrc
```
then add the following two lines to your bashrc:
```
export asfem=/path/to/your/asfem/
export PATH=$PATH:$asfem/bin
```
then save it and open a new terminal, you can run AsFem as follow:
```
asfem -i your_input.i
```
or in the parallel way:
```
mpirun -np 4 asfem -i your_input.i
```

Enjoy it!
