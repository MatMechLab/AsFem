# Basic components

Before we begin the installation, some components are needed, for example:

* gcc g++ gfortran
* cmake
* git
* python3

<span style="color:red">In the following installation, please **do not** use `sudo` or `root`!!! It's only supposed to be used for `apt install`</span>.

**It is highly recommended to install openmpi and PETSc from the source code!**

## Install gcc
You can use the steps below to install a modern GCC compiler if you don't already have one. **Otherwise, go straight to the MPI setup**.

<span style="color:red">**Before you start, you'll need a GCC (the old one) compiler to compile another GCC, therefore if your system doesn't have one, please run(for Ubuntu)**</span>:
```
sudo apt install gcc g++ gfortran cmake git python3 build-essential
```
For ubuntu 20.04, `apt install` will bring you the `GCC9.3.0` package which is **new** and **good** enough for our installation, so you can skip the GCC installation part list below. Anyway, you can also install a newer version, i.e., GCC10.3.0 or even newer.

To begin, you should download and unzip the gcc source code:
```
curl -L -O http://mirrors.concertpass.com/gcc/releases/gcc-10.3.0/gcc-10.3.0.tar.gz

tar -xf gcc-10.3.0.tar.gz
```
or directly go to the official website [gcc](https://gcc.gnu.org/).

Next, we'll need to get the GCC pre-request:
```
cd gcc-10.3.0
./contrib/download_prerequisites
```
we can then config, build and install GCC:
```
./configure --prefix=**your-gcc-install-path** \
--disable-multilib \
--enable-languages=c,c++,fortran,jit \
--enable-checking=release \
--enable-host-shared \
--with-pic
```

```
make -j4
```

```
make install
```

For the following installation, we must configure our bash environment to use our modern GCC:
```
export gcc=your-gcc-instal-path
export PATH=$gcc/bin:$PATH
export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/10.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/10.3.0:$LD_LIBRARY_PATH
```
you can either paste the above settings to your `~/.bashrc` file or put it into a new file, i.e., `~/.asfem-profile`.


## Install mpi
The MPI compiler is required to use the PETSc package for parallelization. **Please skip this stage if your machine already has the MPI compiler installed**.
```
curl -L -O https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz
```
Alternatively, you can get the source code from [openmpi](https://www.open-mpi.org/).

Then you can config, build and install your openmpi (<span style="color:red">if you installed the GCC by compiling the source code rather than using `apt install`, then please do `source ~/.asfem-profile` to set up your environment before running the following step</span>):
```
tar -xf openmpi-4.1.0.tar.gz
cd openmpi-4.1.0
./configure --prefix=*your-path-to-opemmpi*
make -j8
make install
```
After that, you should put the related settings into your bash environment as follows:
```
export MPI_DIR=your-path-to-openmpi

export PATH=$PATH:$MPI_DIR/bin

# please comment out the following four lines until the PETSc has installed !!!
#export CC=mpicc
#export CXX=mpicxx
#export FC=mpif90
#export F90=mpif90

export OMP_NUM_THREADS=1

export C_INCLUDE_PATH=$MPI_DIR/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$MPI_DIR/include:$CPLUS_INCLUDE_PATH
export FPATH=$MPI_DIR/include:$FPATH
export MANPATH=$MPI_DIR/share/man:$MANPATH
export LD_LIBRARY_PATH=$MPI_DIR/lib:$LD_LIBRARY_PATH
```

Now execute the `source ~/.asfem-profile`, and then run:
```
mpirun -np 4 echo "Hi"
```
you should see 4x"Hi" in your terminal.

# Install PETSc
**Before we begin, make sure your GCC and MPI compilers are up to date (by doing `source ~/.asfem-profile`)**, (<span style="color:red">if both of your GCC and MPI are installed from `apt install`, then the newly opened terminal will find the gcc/mpicc compiler, which means you **do not** need to run `source ~/.bashrc` or `source ~/.asfem-profile` !</span>):
```
source ~/.asfem-profile
gcc --version
mpicxx --version
```
Plase keep in mind, these four lines must be comment out in your `~/.asfem-profile`
```
#export CC=mpicc
#export CXX=mpicxx
#export FC=mpif90
#export F90=mpif90
```
We will tell PETSc the path to the openmpi, then it will find the correct compiler.

The PETSc package can be downloaded via(you can change the version number to whatever you like, for example, `petsc-3.16.0.tar.gz`, `petsc-3.13.3.tar.gz`):
```
curl -L -O  https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.15.3.tar.gz
```
or you can download it from the [PETSc website](https://www.mcs.anl.gov/petsc/download/index.html).

For the configuration, one can use:
```
./configure \
--prefix=***your-PETSc-install-path*** \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=1 \
--with-mpi-dir=***your-MPI-install-path*** \
--with-shared-libraries=1 \
--with-cxx-dialect=C++14 \
--with-fortran-bindings=1 \
--with-sowing=0 \
--download-hypre=1 \
--download-fblaslapack=1 \
--download-ptscotch=1 \
--download-metis=1 \
--download-parmetis=1 \
--download-superlu_dist=1 \
--download-scalapack=1 \
--download-mumps=1 \
--download-slepc=1 \
--download-fftw=1 \
--download-hwloc=1 \
COPTFLAGS='-fPIC -fopenmp -O3 ' \
CXXOPTFLAGS='-fPIC -fopenmp -O3 ' \
FOPTFLAGS='-fPIC -fopenmp -O3 ' \
PETSC_DIR=`pwd`
```
you can also use the simplified version (with much fewer external packages) as follows:
```
./configure \
--prefix=***your-PETSc-install-path*** \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=1 \
--with-mpi-dir=***your-MPI-install-path*** \
--with-shared-libraries=1 \
--with-cxx-dialect=C++14 \
--with-fortran-bindings=1 \
--with-sowing=0 \
--download-fblaslapack=1 \
COPTFLAGS='-fPIC -fopenmp -O3 ' \
CXXOPTFLAGS='-fPIC -fopenmp -O3 ' \
FOPTFLAGS='-fPIC -fopenmp -O3 ' \
PETSC_DIR=`pwd`
```
other external packages can be installed via `-download-XXX=1` . Then, by running the `make xxxxx -j8` and `make xxxx install` command, one can finish the installation. Here `xxx` represents the command line which is shown by PETSc in your terminal!

Afterwards, put the related settings into your `~/.asfem-profile` (now, you need to uncomment `export CC=xxx`) as follows:
```
export gcc=***your-path-to-gcc-install-path***
export PETSC_DIR=***your-path-to-petsc-install-dir***
export MPI_DIR=***your-path-to-openmpi***

export PATH=$gcc/bin:$PATH
export PATH=$MPI_DIR/bin:$PATH

export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/10.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/10.3.0:$LD_LIBRARY_PATH

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90

export C_INCLUDE_PATH=$MPI_DIR/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$MPI_DIR/include:$CPLUS_INCLUDE_PATH
export FPATH=$MPI_DIR/include:$FPATH
export MANPATH=$MPI_DIR/share/man:$MANPATH
export LD_LIBRARY_PATH=$MPI_DIR/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
```

Once again, if you installed your GCC and openmpi from `apt install`, then you only need:
```
export PETSC_DIR=***your-path-to-petsc-install-dir***
export MPI_DIR=***your-path-to-openmpi***

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90

export OMP_NUM_THREADS=1
```

# Install AsFem

Download AsFem:
```
git clone https://github.com/yangbai90/AsFem.git
```
if you want to have a try for the `devel` version(unstable), one can use:
```
git clone -b devel https://github.com/yangbai90/AsFem.git
```
then execute `source ~/.asfem-profile`, and then we can get the makefile via:
```
cmake CMakeLists.txt
```
afterwards, we can make the `asfem` by executing:
```
make -j4
```
the executable file `asfem` can be found in the `AsFem/bin` folder.
