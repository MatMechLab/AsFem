# Basic components

Before we begin the installation, some components are needed, for example:

* gcc
* cmake
* git



## Install gcc
You can use the steps below to install a modern GCC compiler if you don't already have one. **Otherwise, go straight to the MPI setup**.

To begin, you should download and unzip the gcc source code:
```
curl -L -O http://mirrors.concertpass.com/gcc/releases/gcc-9.3.0/gcc-9.3.0.tar.gz

tar -xf gcc-9.3.0.tar.gz
```
or directly go to the official website [gcc](https://gcc.gnu.org/).

Next, we'll need to get the GCC pre-request:
```
cd gcc-9.3.0
./contrib/download_prerequisites
```
we can then config, build and install GCC:
```
./configure --prefix=your-gcc-install-path \
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

Following the installation, we must configure our bash environment to use our modern GCC:
```
export gcc=your-gcc-instal-path
export PATH=$gcc/bin:$PATH
export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/9.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/9.3.0:$LD_LIBRARY_PATH
```
you can either paste the above settings to your `~/.bashrc` file or put it into a new file, i.e. `~/.asfem-profile`.


## Install mpi
The MPI compiler is required to use the PETSc package for parallelization. **Please skip this stage if your machine already has the MPI compiler installed**.
```
curl -L -O https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.0.tar.gz
```
Alternatively, you can get the source code from [openmpi](https://www.open-mpi.org/).

Then you can config, build and install your openmpi:
```
tar -xf openmpi-4.1.0.tar.gz
cd openmpi-4.1.0
./configure --prefix=*your-path-to-opemmpi*
make -j8
make install
```
After that, you should put the related settings in your bash environment as follows:
```
export mpi=your-path-to-openmpi

export PATH=$PATH:$mpi/bin

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export OMP_NUM_THREADS=1

export C_INCLUDE_PATH=$mpi/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$mpi/include:$CPLUS_INCLUDE_PATH
export FPATH=$mpi/include:$FPATH
export MANPATH=$mpi/share/man:$MANPATH
export LD_LIBRARY_PATH=$mpi/lib:$LD_LIBRARY_PATH
```

# Install PETSc
**Before we begin, make sure your GCC and MPI compilers are up to date**:
```
source ~/.asfem-profile
gcc --version
mpicxx --version
```

The PETSc package can be downloaded via:
```
curl -L -O  https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.14.3.tar.gz
```
or you can download it from the [PETSc website](https://www.mcs.anl.gov/petsc/download/index.html).

For the configuration, one can use:
```
./configure \
--prefix=your-PETSc-install-path \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=1 \
--with-mpi-dir=your-MPI-install-path \
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
other external packages can be installed via `-download-XXX=1` . Then, by running the `make xxxxx -j8` and `make xxxx install` command, one can finish the installation. Here `xxx` represents the command line which is shown by PETSc in your terminal!

# Install AsFem

Download AsFem:
```
git clone https://github.com/yangbai90/AsFem.git
```
if you want the have a try for the `devel` version(unstable), one can use:
```
git clone -b devel https://github.com/yangbai90/AsFem.git
```
then we can make `asfem` via:
```
make -j4
```
the executable file `asfem` can be found in the `AsFem/bin` folder.
