cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 17)

###############################################
### Set your Eigen/PETSc/MPI path here      ###
###############################################
set(PETSC_DIR "/home/by/Programs/petsc/3.12.5")
set(MPI_DIR "/home/by/Programs/mpich/3.3.2")
set(Eigen "/home/by/Programs/eigen")

###############################################
### For your PETSc library file  (edit)     ###
###############################################
include_directories("${PETSC_DIR}/include")
link_libraries("${PETSC_DIR}/lib/libpetsc.so")
###############################################
### For your MPI library file (edit)        ###
###############################################
include_directories("${MPI_DIR}/include")
link_libraries("${MPI_DIR}/lib/libmpi.so")

# For Eigen
include_directories("${Eigen}")

###############################################
### For debug mode                          ###
###############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O2 -fopenmp")

###############################################
### For release mode                        ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG -Wno-deprecated")
# if you have mumps installed, please use the following line
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG -DHASMUMPS")


###############################################
### For release mode in Windows (VS2019)    ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W1 /arch:AVX")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /openmp")


###############################################
### Do not edit the following two lines !!! ###
###############################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#############################################################
#############################################################
### For beginners, please don't edit the following line!  ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
#############################################################
#############################################################
# For Welcome header file and main.cpp
set(inc include/Welcome.h)
set(src src/main.cpp)

#############################################################
### For String utils                                      ###
#############################################################
set(inc ${inc} include/Utils/StringUtils.h)
set(src ${src} src/Utils/StringUtils.cpp)

#############################################################
### For message printer utils                             ###
#############################################################
set(inc ${inc} include/Utils/MessagePrinter.h)
set(src ${src} src/Utils/MessagePrinter.cpp)

#############################################################
### For inputystem                                        ###
#############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp src/InputSystem/ReadInputFile.cpp)
set(src ${src} src/InputSystem/ReadMeshBlock.cpp)
set(src ${src} src/InputSystem/ReadDofsBlock.cpp)
set(src ${src} src/InputSystem/ReadElmtBlock.cpp)
set(src ${src} src/InputSystem/ReadMateBlock.cpp)

#############################################################
### For Mesh class                                        ###
#############################################################
### For mesh io class
set(inc ${inc} include/Mesh/MeshIOBase.h include/Mesh/Gmsh2IO.h)
set(src ${src} src/Mesh/MeshIO/Gmsh2IOReadMesh.cpp src/Mesh/MeshIO/Gmsh2IOGetFuns.cpp)
set(src ${src} include/Mesh/Gmsh4IO.h)
set(src ${src} src/Mesh/MeshIO/Gmsh4IOReadMesh.cpp src/Mesh/MeshIO/Gmsh4IOGetFuns.cpp)
set(src ${src} include/Mesh/MeshIO.h src/Mesh/MeshIO/MeshIO.cpp src/Mesh/MeshIO/MeshIOJudgeFun.cpp)
### define the basic mesh type and geometry
set(inc ${inc} include/Mesh/MeshType.h include/Mesh/Nodes.h)
### for the bulk mesh class
set(inc ${inc} include/Mesh/LagrangeMesh.h)
set(src ${src} src/Mesh/LagrangeMesh.cpp src/Mesh/CreateLagrangeMesh.cpp)
set(src ${src} src/Mesh/Create1DLagrangeMesh.cpp)
set(src ${src} src/Mesh/Create2DLagrangeMesh.cpp)
set(src ${src} src/Mesh/Create3DLagrangeMesh.cpp)
set(src ${src} src/Mesh/SaveLagrangeMesh.cpp)
set(src ${src} src/Mesh/LagrangeMeshPrintInfo.cpp)

### for the final mesh class
set(inc ${inc} include/Mesh/Mesh.h)
set(src ${src} src/Mesh/Mesh.cpp)


#############################################################
### For DofHandler class                                  ###
#############################################################
set(inc ${inc} include/DofHandler/BulkDofHandler.h)
set(src ${src} src/DofHandler/BulkDofHandler.cpp)
set(src ${src} src/DofHandler/BulkDofHandlerFuns.cpp)
set(src ${src} src/DofHandler/BulkDofHandlerGetFuns.cpp)
### for final dof class
set(inc ${inc} include/DofHandler/DofHandler.h)
set(src ${src} src/DofHandler/DofHandler.cpp)

#############################################################
### For Element system in AsFem                         ###
#############################################################
set(inc ${inc} include/ElmtSystem/ElmtType.h)
set(inc ${inc} include/ElmtSystem/ElmtBlock.h)
set(inc ${inc} include/ElmtSystem/ElmtSystem.h)
set(src ${src} src/ElmtSystem/ElmtSystem.cpp)

#############################################################
### For Materials system in AsFem                         ###
#############################################################
set(inc ${inc} include/MateSystem/MateType.h)
set(inc ${inc} include/MateSystem/MateBlock.h)
set(inc ${inc} include/MateSystem/MateSystem.h)
set(src ${src} src/MateSystem/MateSystem.cpp)

##################################################
add_executable(asfem ${inc} ${src})

