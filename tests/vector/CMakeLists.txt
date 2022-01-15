cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 17)

if(UNIX)
    message ("We are running on linux system ...")
elseif(MSVC)
    message("We are running on windows system (MSVC) ...")
endif()

###############################################
### Set your PETSc/MPI path here or bashrc  ###
### The only things to modify is the        ###
### following two lines(PETSC/MPI_DIR)      ###
###############################################


if(EXISTS $ENV{MPI_DIR})
    set(MPI_DIR $ENV{MPI_DIR})
    message("MPI dir is: ${MPI_DIR}")
else()
    message (WARNING "MPI location (MPI_DIR) is not defined in your PATH, AsFem will use the one defined in CMakeLists.txt")
    set(MPI_DIR "/home/by/Programs/openmpi/4.1.0")
    message("MPI dir set to be: ${MPI_DIR}")
    message (WARNING "If the path is not correct, you should modify line-24 in your CMakeLists.txt")
endif()


if(EXISTS $ENV{PETSC_DIR})
    set(PETSC_DIR $ENV{PETSC_DIR})
    message("PETSC dir is: ${PETSC_DIR}")
else()
    message (WARNING "PETSc location (PETSC_DIR) is not defined in your PATH, AsFem will use the one defined in CMakeLists.txt")
    set(PETSC_DIR "/home/by/Programs/petsc/3.14.3")
    message("PETSc dir set to be:${PETSC_DIR}")
    message (WARNING "If the path is not correct, you should modify line-35 in your CMakeLists.txt")
endif()

get_filename_component(ASFEM_DIR ../../ ABSOLUTE)
message("AsFem dir is:${ASFEM_DIR}")

###############################################
### For include files of PETSc and mpi      ###
###############################################
include_directories("${PETSC_DIR}/include")
include_directories("${MPI_DIR}/include")
if(UNIX)
    link_libraries("${PETSC_DIR}/lib/libpetsc.so")
    link_libraries("${MPI_DIR}/lib/libmpi.so")
elseif(MSVC)
    link_libraries("${PETSC_DIR}/lib/libpetsc.lib")
endif()

###############################################
# For Eigen                                 ###
###############################################
include_directories("${ASFEM_DIR}/external/eigen")


###############################################
### set debug or release mode               ###
###############################################
if (CMAKE_BUILD_TYPE STREQUAL "")
    # user should use -DCMAKE_BUILD_TYPE=Release[Debug] option
    set (CMAKE_BUILD_TYPE "Debug")
endif ()

###############################################
### For linux platform                      ###
###############################################
if(UNIX)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O2 -fopenmp")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG")
    else()
        message (FATAL_ERROR "Unknown compiler flags (CMAKE_CXX_FLAGS)")
    endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W1 /arch:AVX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /openmp")
endif()

message("AsFem will be compiled in ${CMAKE_BUILD_TYPE} mode !")


###############################################
### Do not edit the following two lines !!! ###
###############################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(${ASFEM_DIR}/include)

#############################################################
#############################################################
### For beginners, please don't edit the following line!  ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
#############################################################
#############################################################
# For Welcome header file and main.cpp
set(inc "")
set(src test.cpp)


#############################################################
### For message printer utils                             ###
#############################################################
set(inc ${inc} ${ASFEM_DIR}/include/Utils/MessagePrinter.h ${ASFEM_DIR}/include/Utils/MessageColor.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MessagePrinter.cpp)
#############################################################
### For mathematic utils (vector and tensors, etc...)     ###
#############################################################
set(inc ${inc} ${ASFEM_DIR}/include/Utils/Vector3d.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MathUtils/Vector3d.cpp)
### for rank-2 tensor
set(inc ${inc} ${ASFEM_DIR}/include/Utils/RankTwoTensor.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MathUtils/RankTwoTensor.cpp)
### for rank-4 tensor
set(inc ${inc} ${ASFEM_DIR}/include/Utils/RankFourTensor.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MathUtils/RankFourTensor.cpp)
### for MatrixXd and VectorXd
set(inc ${inc} ${ASFEM_DIR}/include/Utils/VectorXd.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MathUtils/VectorXd.cpp)
set(inc ${inc} ${ASFEM_DIR}/include/Utils/MatrixXd.h)
set(src ${src} ${ASFEM_DIR}/src/Utils/MathUtils/MatrixXd.cpp)
### for general mathematic functions
set(inc ${inc} ${ASFEM_DIR}/include/Utils/MathFuns.h)


##################################################
add_executable(asfem-test ${inc} ${src})


##################################################
### Following lines are used by vim            ###
### you can delete all of them                 ###
##################################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PETSC_DIR}/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_DIR}/include")

