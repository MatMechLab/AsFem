# System requirement

* Linux
* Windows (only works in Cygwin/Mingw, <span style="color:red">**without**</span> **VisualStudio**)
* MacOX (it should work, but we haven't tested it yet!)


# Basic components

Before we start the installation, some components are required, for example:

* gcc g++ gfortran
* cmake
* git
* python3

<span style="color:red">In the following installation, please **do not** use `sudo` or `root`!!! It's only supposed to be used for system installation like `sudo apt install`, `sudo zypper install`, `sudo yum install`, and so on</span>.

**It is highly recommended to install the openmpi and PETSc packages from the source code!**

## Install gcc
One can follow the steps list below to install a modern GCC compiler if one don't already have one. **Otherwise, please go straight to the MPI setup**.

<span style="color:red">**Before you start, you'll need a GCC (the old one) compiler to compile another GCC, therefore if your system doesn't have one, please run(for Ubuntu)**</span>:
```
sudo apt install gcc g++ gfortran cmake git python3 build-essential
```
For ubuntu 22.04, `apt install` will bring you the `GCC11.2.0` package which is **new** and **good** enough for our installation, so you can skip the GCC installation part list below. Anyway, you can also install a newer version, i.e., GCC12.1.0 or even newer.

<span style="color:red">**The installation itself is not limited to Ubuntu system, for other linux distributions, one need to change the command to the consistent one. For instance, `sudo zypper`, `sudo yum`, `sudo dnf`...**</span>:

To begin, one should download and unzip the gcc source code:
```
curl -L -O http://mirrors.concertpass.com/gcc/releases/gcc-11.3.0/gcc-11.3.0.tar.gz

tar -xf gcc-11.3.0.tar.gz
```
or directly go to the official website [gcc](https://gcc.gnu.org/).

Next, we'll need to download the GCC pre-request as follows:
```
cd gcc-11.3.0
./contrib/download_prerequisites
```
we can then configure, build and install GCC:
```
./configure --prefix=**your-gcc-install-path** \
--disable-multilib \
--enable-languages=c,c++,fortran,jit \
--enable-checking=release \
--enable-host-shared \
--with-pic
```
where `**your-gcc-install-path**` represents the installation path on **your own computer**.

Then one can execute:
```
make -j4
```
and
```
make install
```

In the following installation, we must configure our bash environment for using our modern GCC(if you have installed your gcc via `sudo apt` or `sudo yum` or `sudo zypper` or whatever `sudo x`, please skip this step, because your GCC is ready!):
```
export gcc=your-gcc-instal-path
export PATH=$gcc/bin:$PATH
export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/11.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/11.3.0:$LD_LIBRARY_PATH
```
you can either paste the above settings to your `~/.bashrc` file or put it into a new file, i.e., `~/.asfem-profile`.


## Install mpi
The MPI compiler is required to use the PETSc package for parallelization. **Please skip this stage if your machine already has the MPI compiler installed**.
```
curl -L -O https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz
```
Alternatively, you can get the source code from [openmpi](https://www.open-mpi.org/).

Then you can config, build and install your openmpi (<span style="color:red">if you installed the GCC by compiling the source code rather than using `sudo apt install`, then please do `source ~/.asfem-profile` to set up your environment before running the following step</span>):
```
tar -xf openmpi-4.1.4.tar.gz
cd openmpi-4.1.4
./configure --prefix=*your-path-to-opemmpi*
make -j8
make install
```
again, `*your-path-to-opemmpi*` should be your own installation path.

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

Now one can execute the `source ~/.asfem-profile` command, and then do:
```
mpirun -np 4 echo "Hi"
```
you should see 4x"Hi" in your terminal.

# Install PETSc
**Before we begin, make sure your GCC and MPI compilers are up to date (by doing `source ~/.asfem-profile`)**, (<span style="color:red">if both of your GCC and MPI are installed from `sudo`, then the newly opened terminal will find the gcc/mpicc compiler, which means you **do not** need to execute `source ~/.bashrc` or `source ~/.asfem-profile` !</span>):
```
source ~/.asfem-profile
gcc --version
mpicxx --version
```
Plase keep in mind, these four lines must be **uncommented** in your `~/.asfem-profile`
```
#export CC=mpicc
#export CXX=mpicxx
#export FC=mpif90
#export F90=mpif90
```
We will tell PETSc the path to the openmpi, then it will find the correct compiler.

The PETSc package can be downloaded via (you can change the version number to whatever you like, for example, `petsc-3.17.2.tar.gz`, `petsc-3.16.3.tar.gz`):
```
curl -L -O  https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.17.2.tar.gz
```
or you can download it from the [PETSc website](https://www.mcs.anl.gov/petsc/download/index.html).

For the configuration, one can use:
```
./configure \
--prefix=***your-PETSc-install-path*** \
--with-debugging=0 \
--with-ssl=0 \
--with-pic=1 \
--with-openmp=0 \
--with-mpi-dir=***your-MPI-install-path*** \
--with-shared-libraries=1 \
--with-cxx-dialect=C++14 \
--with-fortran-bindings=1 \
--with-sowing=0 \
--download-fblaslapack=1 \
COPTFLAGS='-fPIC -O3 -march=native -mtune=native ' \
CXXOPTFLAGS='-fPIC -O3 -march=native -mtune=native ' \
FOPTFLAGS='-fPIC -O3 -march=native -mtune=native ' \
PETSC_DIR=`pwd`
```
once again, `***your-PETSc-install-path***` should be your PETSc installation path, and `***your-MPI-install-path***` is the path of your MPI. Then, by running the `make xxxx -j8` and `make xxxx install` command, one can finish the installation. Here `xxx` represents the command line which is shown by PETSc in your terminal!

If `-march=native` and `-mtune=native ` flag doesn't work on your system, then please remove it.

Afterwards, put the related settings into your `~/.asfem-profile` (now, you need to uncomment `export CC=xxx`) as follows:
```
export gcc=***your-path-to-gcc-install-path***
export PETSC_DIR=***your-path-to-petsc-install-dir***
export MPI_DIR=***your-path-to-openmpi***

export PATH=$gcc/bin:$PATH
export PATH=$MPI_DIR/bin:$PATH

export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/11.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/11.3.0:$LD_LIBRARY_PATH

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

Once again, if you installed your GCC and openmpi from `sudo`, then you only need:
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
git clone https://github.com/M3Group/AsFem.git
```
if one wants to have a try for the `devel` version (**unstable**), one can use:
```
git clone -b devel https://github.com/M3Group/AsFem.git
```
then execute `source ~/.asfem-profile`, and then we can get the makefile via:
```
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
```
afterwards, we can make the `asfem` by executing:
```
make -j4
```
the executable file `asfem` can be found in the `AsFem/bin` folder.

It is also highly recommended to put `asfem` into your `PATH` as follows:
```
export asfem=**your-path-to-asfem**
export PATH=$PATH:$asfem/bin
```
then one can easily run the AsFem job from the terminal as follows:
```
asfem -i yourjob.i
mpirun -np 16 asfem -i yourjob.i
```

That's all, enjoy!
