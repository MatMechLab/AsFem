# System requirement

* Linux
* Windows (only works in Cygwin/Mingw, <span style="color:red">**without**</span> **VisualStudio**. <span style="color:red">**It is not recommended!!!**</span>)
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

If you don't already have a modern GCC compiler, you can use the steps listed below to install one. **Otherwise, please proceed directly to the MPI setup.**

<span style="color:red">**To compile a new GCC compiler, you will first need an old GCC compiler. If your system does not have one, please execute the following command (for Ubuntu) before proceeding.**</span>:
```
sudo apt install gcc g++ gfortran cmake git python3 build-essential
```
If you're using Ubuntu 22.04, the apt install command will provide you with the `GCC11.3.0` package, which is a suitable and up-to-date choice for our installation. As a result, you can skip the list of steps for installing GCC below. However, if you prefer, you can also install a more recent version of GCC, such as `GCC12.2.0` or even later.

<span style="color:red">**The installation process is not exclusive to Ubuntu systems. For other Linux distributions, you may need to adjust the command to match the appropriate package manager. For example, you may need to use `sudo zypper`, `sudo yum`, `sudo dnf`, and so on.**</span>:

To get started, you will need to download and extract the source code for GCC:
```
curl -L -O http://mirrors.concertpass.com/gcc/releases/gcc-11.3.0/gcc-11.3.0.tar.gz

tar -xf gcc-11.3.0.tar.gz
```
Alternatively, you can also directly visit the official [GCC](https://gcc.gnu.org/) website to obtain the source code.

Next, you will need to download the GCC prerequisites by following these steps:
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
where <span style="color:red">`**your-gcc-install-path**`</span> represents the installation path on **your own computer**.

Then one can execute:
```
make -j4
```
and
```
make install
```

During the installation process, it is necessary to configure your bash environment to use the modern GCC compiler. However, if you have installed GCC using `sudo apt`, `sudo yum`, `sudo zypper`, or any similar command, you can skip this step since your GCC is already set up and ready to use.
```
export gcc=your-gcc-instal-path
export PATH=$gcc/bin:$PATH
export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/11.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/11.3.0:$LD_LIBRARY_PATH
```
You can add the aforementioned settings to your `~/.bashrc` file, or create a new file named `~/.asfem-profile` and add them there.


## Install mpi

**If your machine already has an MPI compiler installed, you can skip this step**. Otherwise, you will need to install it in order to use the PETSc package for parallelization.

You can download the OpenMPI source code using the following command:
```
curl -L -O https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz
```
Alternatively, you can get the source code from [openmpi](https://www.open-mpi.org/).

Next, you can configure, build, and install OpenMPI. <span style="color:red">If you installed GCC by compiling the source code rather than using `sudo apt install`, please execute `source ~/.asfem-profile` to set up your environment before proceeding with the following step</span>.
```
tar -xf openmpi-4.1.5.tar.gz
cd openmpi-4.1.5
./configure --prefix=*your-path-to-opemmpi*
make -j8
make install
```
once again, `*your-path-to-opemmpi*` should be <span style="color:red">your own installation path</span>.

Once you have installed OpenMPI, you will need to add the necessary settings to your bash environment as follows:
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
**Before proceeding, ensure that your GCC and MPI compilers are up to date by executing `source ~/.asfem-profile`**, (<span style="color:red">If you have installed both GCC and MPI using sudo commands, then the newly opened terminal will automatically recognize the gcc and mpicc compilers, and you will not need to execute `source ~/.bashrc` or `source ~/.asfem-profile` !</span>):
```
source ~/.asfem-profile
gcc --version
mpicxx --version
```
Plase keep in mind, these four lines should be **commented** in your `~/.asfem-profile` (<span style="color:red">You can also **uncomment** these lines, but please make sure that the `MPI_DIR` is set correctly to match the selected mpicxx compiler.</span>)
```
#export CC=mpicc
#export CXX=mpicxx
#export FC=mpif90
#export F90=mpif90
```
By specifying the path to OpenMPI (mpicxx,mpicc, etc.), PETSc will be able to locate the appropriate compiler.

You can download the PETSc package using the following command (you can change the version number to any desired version, such as `petsc-3.18.2.tar.gz` or `petsc-3.19.1.tar.gz`):
```
curl -L -O  https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.18.6.tar.gz
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
Once again, please note that `***your-PETSc-install-path***` should be replaced with the actual path ***where you want to install PETSc***, and `***your-MPI-install-path***` should be replaced with the path to your MPI installation. <span style="color:red">Please ensure that the `--with-mpi-dir` option is set up correctly</span>.

Finally, to complete the installation, run the command shown by PETSc in your terminal using the `make xxxx -j8` and `make xxxx install` commands. Replace ***`xxx`*** with the appropriate command line shown in your terminal(<span style="color:red">Note that the commands for your particular installation **may be different**, so please **do not copy and modify the commands shown here**. Instead, copy the commands specific to your own installation **from your terminal**!</span>)

If the `-march=native` and `-mtune=native` flags do not work on your system, you should remove them from the compiler flags.

It is important to note that if you wish to use the **`lu`** preconditioner or the **direct solver**, you will need to add one of the following lines for the lu solver:
```
--download-superlu_dist=1
--download-mumps=1
```

Afterwards, add the necessary settings to your `~/.asfem-profile`. Please note that you will need to uncomment `export CC=xxx` options. Here is an example of a complete `~/.asfem-profile` file that includes the necessary settings for using PETSc.
```
export gcc=***your-path-to-gcc-install-path***
export PETSC_DIR=***your-path-to-petsc-install-dir***
export MPI_DIR=***your-path-to-openmpi***

export PATH=$gcc/bin:$PATH
export PATH=$MPI_DIR/bin:$PATH

export LD_LIBRARY_PATH=$gcc/lib64:$gcc/lib:$gcc/lib/gcc/x86_64-pc-linux-gnu/11.3.0:$gcc/libexec/gcc/x86_64-pc-linux-gnu/11.3.0:$LD_LIBRARY_PATH

export CC=$MPI_DIR/bin/mpicc
export CXX=$MPI_DIR/bin/mpicxx
export FC=$MPI_DIR/bin/mpif90
export F90=$MPI_DIR/bin/mpif90

export C_INCLUDE_PATH=$MPI_DIR/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$MPI_DIR/include:$CPLUS_INCLUDE_PATH
export FPATH=$MPI_DIR/include:$FPATH
export MANPATH=$MPI_DIR/share/man:$MANPATH
export LD_LIBRARY_PATH=$MPI_DIR/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
```

Once again, if you have installed GCC and OpenMPI using `sudo` commands, you only need to add the PETSc settings to your `~/.asfem-profile`. Here are some sample commands for configuring your profile:
```
export PETSC_DIR=***your-path-to-petsc-install-dir***
export MPI_DIR=***your-path-to-openmpi***

export CC=$MPI_DIR/bin/mpicc
export CXX=$MPI_DIR/bin/mpicxx
export FC=$MPI_DIR/bin/mpif90
export F90=$MPI_DIR/bin/mpif90

export OMP_NUM_THREADS=1
```

# Install AsFem

Download AsFem:
```
git clone https://github.com/M3Group/AsFem.git
```
if you wish to try the devel version (which is considered to be **unstable**), you can use the following command:
```
git clone -b devel https://github.com/M3Group/AsFem.git
```

After executing source ~/.asfem-profile, you can obtain the `Makefile` using the following command:
```
cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release
```
After obtaining the `Makefile`, you can build `asfem` by executing the following command:
```
make -j4
```
Once you have built `asfem`, you can find the executable file in the `AsFem/bin` folder.

It is highly recommended to add the `AsFem/bin` directory to your `PATH` environment variable as follows:
```
export asfem=**your-path-to-asfem**
export PATH=$PATH:$asfem/bin
```
After adding the `AsFem/bin` directory to your PATH environment variable, you can easily run an AsFem job from the terminal using the following command:
```
asfem -i yourjob.json
mpirun -np 16 asfem -i yourjob.json
```

That's all, enjoy!
