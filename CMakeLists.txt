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
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")


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
set(inc ${inc} include/Utils/MessagePrinter.h include/Utils/MessageColor.h)
set(src ${src} src/Utils/MessagePrinter.cpp)

#############################################################
### For mathematic utils (vector and tensors, etc...)     ###
#############################################################
set(inc ${inc} include/Utils/Vector3d.h)
set(src ${src} src/Utils/MathUtils/Vector3d.cpp)
### for rank-2 tensor
set(inc ${inc} include/Utils/RankTwoTensor.h)
set(src ${src} src/Utils/MathUtils/RankTwoTensor.cpp)
### for rank-4 tensor
set(inc ${inc} include/Utils/RankFourTensor.h)
set(src ${src} src/Utils/MathUtils/RankFourTensor.cpp)
### for MatrixXd and VectorXd
set(inc ${inc} include/Utils/VectorXd.h)
set(src ${src} src/Utils/MathUtils/VectorXd.cpp)
set(inc ${inc} include/Utils/MatrixXd.h)
set(src ${src} src/Utils/MathUtils/MatrixXd.cpp)
### for general mathematic functions
set(inc ${inc} include/Utils/MathFuns.h)

#############################################################
### For inputystem                                        ###
#############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp src/InputSystem/ReadInputFile.cpp)
set(inc ${inc} include/InputSystem/SingleBlockReader.h)
### for mesh block
set(inc ${inc} include/InputSystem/MeshBlockReader.h)
set(src ${src} src/InputSystem/MeshBlockReader.cpp)
### for dof block
set(inc ${inc} include/InputSystem/DofsBlockReader.h)
set(src ${src} src/InputSystem/DofsBlockReader.cpp)
### for elmts block
set(inc ${inc} include/InputSystem/ElmtsBlockReader.h)
set(src ${src} src/InputSystem/ElmtsBlockReader.cpp)
### for mate block
set(inc ${inc} include/InputSystem/MatesBlockReader.h)
set(src ${src} src/InputSystem/MatesBlockReader.cpp)
### for bcs block
set(inc ${inc} include/InputSystem/BCsBlockReader.h)
set(src ${src} src/InputSystem/BCsBlockReader.cpp)
### for ics block
set(inc ${inc} include/InputSystem/ICsBlockReader.h)
set(src ${src} src/InputSystem/ICsBlockReader.cpp)
### for qpoint block
set(inc ${inc} include/InputSystem/QPointBlockReader.h)
set(src ${src} src/InputSystem/QPointBlockReader.cpp)
### for output block
set(inc ${inc} include/InputSystem/OutputBlockReader.h)
set(src ${src} src/InputSystem/OutputBlockReader.cpp)
### for postprocess block 
set(inc ${inc} include/InputSystem/PostprocessBlockReader.h)
set(src ${src} src/InputSystem/PostprocessBlockReader.cpp)
### for projection block
set(inc ${inc} include/InputSystem/ProjectionBlockReader.h)
set(src ${src} src/InputSystem/ProjectionBlockReader.cpp)
### for nonlinear solver block
set(inc ${inc} include/InputSystem/NonlinearSolverBlockReader.h)
set(src ${src} src/InputSystem/NonlinearSolverBlockReader.cpp)
### for time stepping block
set(inc ${inc} include/InputSystem/TimeSteppingBlockReader.h)
set(src ${src} src/InputSystem/TimeSteppingBlockReader.cpp)
### for FEJob block
set(inc ${inc} include/InputSystem/FEJobBlockReader.h)
set(src ${src} src/InputSystem/FEJobBlockReader.cpp)


#############################################################
### For Mesh class                                        ###
#############################################################
### For mesh io class
### for gmsh2 io
set(inc ${inc} include/Mesh/MeshIOBase.h include/Mesh/Gmsh2IO.h)
set(src ${src} src/Mesh/MeshIO/Gmsh2IOReadMesh.cpp src/Mesh/MeshIO/Gmsh2IOGetFuns.cpp)
### for gmsh4 io
set(src ${src} include/Mesh/Gmsh4IO.h)
set(src ${src} src/Mesh/MeshIO/Gmsh4IOReadMesh.cpp src/Mesh/MeshIO/Gmsh4IOGetFuns.cpp)
set(src ${src} src/Mesh/MeshIO/Gmsh4IOUtilFuns.cpp)
set(src ${src} include/Mesh/MeshIO.h src/Mesh/MeshIO/MeshIO.cpp src/Mesh/MeshIO/MeshIOJudgeFun.cpp)
### for abaqus io
set(inc ${inc} include/Mesh/AbaqusIO.h)
set(src ${src} src/Mesh/MeshIO/AbaqusIOGetFuns.cpp)
set(src ${src} src/Mesh/MeshIO/AbaqusIOReadMesh.cpp)
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
set(src ${src} src/DofHandler/BulkDofHandlerCreateMap.cpp)
### for final dof class
set(inc ${inc} include/DofHandler/DofHandler.h)
set(src ${src} src/DofHandler/DofHandler.cpp)

#############################################################
### For Element system in AsFem                           ###
#############################################################
set(inc ${inc} include/ElmtSystem/ElmtType.h)
set(inc ${inc} include/ElmtSystem/ElmtBlock.h)
set(inc ${inc} include/ElmtSystem/ElmtSystem.h)
set(src ${src} src/ElmtSystem/ElmtSystem.cpp)
# for bulk element system
set(inc ${inc} include/ElmtSystem/BulkElmtSystem.h)
set(src ${src} src/ElmtSystem/BulkElmtSystem.cpp)
set(src ${src} src/ElmtSystem/RunBulkElmtLibs.cpp)
### For bulk element base class
set(inc ${inc} include/ElmtSystem/BulkElmtBase.h)
### For the data structure used by local element calc
set(inc ${inc} include/ElmtSystem/LocalElmtData.h)
### For linear poisson element
set(inc ${inc} include/ElmtSystem/PoissonElmt.h)
set(src ${src} src/ElmtSystem/PoissonElmt.cpp)
### For diffusion elmt
set(inc ${inc} include/ElmtSystem/DiffusionElmt.h)
set(src ${src} src/ElmtSystem/DiffusionElmt.cpp)
### For Mechanics element implementation
set(inc ${inc} include/ElmtSystem/MechanicsElmt.h)
set(src ${src} src/ElmtSystem/MechanicsElmt.cpp)
### For CahnHilliard element implementation
set(inc ${inc} include/ElmtSystem/CahnHilliardElmt.h)
set(src ${src} src/ElmtSystem/CahnHilliardElmt.cpp)
### For Miehe's phase field fracture element
set(inc ${inc} include/ElmtSystem/MieheFractureElmt.h)
set(src ${src} src/ElmtSystem/MieheFractureElmt.cpp)
### For Kobayashi's dendrite element
set(inc ${inc} include/ElmtSystem/KobayashiElmt.h)
set(src ${src} src/ElmtSystem/KobayashiElmt.cpp)
### For stress diffusion element
set(inc ${inc} include/ElmtSystem/StressDiffusionElmt.h)
set(src ${src} src/ElmtSystem/StressDiffusionElmt.cpp)
### For diffusion fracture element
set(inc ${inc} include/ElmtSystem/DiffusionFractureElmt.h)
set(src ${src} src/ElmtSystem/DiffusionFractureElmt.cpp)
### For mechanically coupled CahnHilliard element
set(inc ${inc} include/ElmtSystem/MechanicsCahnHilliardElmt.h)
set(src ${src} src/ElmtSystem/MechanicsCahnHilliardElmt.cpp)
### For AllenCahn fracture element
set(inc ${inc} include/ElmtSystem/AllenCahnFractureElmt.h)
set(src ${src} src/ElmtSystem/AllenCahnFractureElmt.cpp)
### For wave element
set(inc ${inc} include/ElmtSystem/WaveElmt.h)
set(src ${src} src/ElmtSystem/WaveElmt.cpp)
### for thermal element
set(inc ${inc} include/ElmtSystem/ThermalElmt.h)
set(src ${src} src/ElmtSystem/ThermalElmt.cpp)
### for UEL1
set(inc ${inc} include/ElmtSystem/User1Elmt.h)
set(src ${src} src/ElmtSystem/User1Elmt.cpp)
### for UEL2
set(inc ${inc} include/ElmtSystem/User2Elmt.h)
set(src ${src} src/ElmtSystem/User2Elmt.cpp)
### for UEL3
set(inc ${inc} include/ElmtSystem/User3Elmt.h)
set(src ${src} src/ElmtSystem/User3Elmt.cpp)
### for UEL4
set(inc ${inc} include/ElmtSystem/User4Elmt.h)
set(src ${src} src/ElmtSystem/User4Elmt.cpp)
### for UEL5
set(inc ${inc} include/ElmtSystem/User5Elmt.h)
set(src ${src} src/ElmtSystem/User5Elmt.cpp)
### for UEL6
set(inc ${inc} include/ElmtSystem/User6Elmt.h)
set(src ${src} src/ElmtSystem/User6Elmt.cpp)
### for UEL7
set(inc ${inc} include/ElmtSystem/User7Elmt.h)
set(src ${src} src/ElmtSystem/User7Elmt.cpp)
### for UEL8
set(inc ${inc} include/ElmtSystem/User8Elmt.h)
set(src ${src} src/ElmtSystem/User8Elmt.cpp)
### for UEL9
set(inc ${inc} include/ElmtSystem/User9Elmt.h)
set(src ${src} src/ElmtSystem/User9Elmt.cpp)
### for UEL10
set(inc ${inc} include/ElmtSystem/User10Elmt.h)
set(src ${src} src/ElmtSystem/User10Elmt.cpp)
### for UEL11
set(inc ${inc} include/ElmtSystem/User11Elmt.h)
set(src ${src} src/ElmtSystem/User11Elmt.cpp)
### for UEL12
set(inc ${inc} include/ElmtSystem/User12Elmt.h)
set(src ${src} src/ElmtSystem/User12Elmt.cpp)
### for UEL13
set(inc ${inc} include/ElmtSystem/User13Elmt.h)
set(src ${src} src/ElmtSystem/User13Elmt.cpp)
### for UEL14
set(inc ${inc} include/ElmtSystem/User14Elmt.h)
set(src ${src} src/ElmtSystem/User14Elmt.cpp)
### for UEL15
set(inc ${inc} include/ElmtSystem/User15Elmt.h)
set(src ${src} src/ElmtSystem/User15Elmt.cpp)
### for UEL16
set(inc ${inc} include/ElmtSystem/User16Elmt.h)
set(src ${src} src/ElmtSystem/User16Elmt.cpp)
### for UEL17
set(inc ${inc} include/ElmtSystem/User17Elmt.h)
set(src ${src} src/ElmtSystem/User17Elmt.cpp)
### for UEL18
set(inc ${inc} include/ElmtSystem/User18Elmt.h)
set(src ${src} src/ElmtSystem/User18Elmt.cpp)
### for UEL19
set(inc ${inc} include/ElmtSystem/User19Elmt.h)
set(src ${src} src/ElmtSystem/User19Elmt.cpp)
### for UEL20
set(inc ${inc} include/ElmtSystem/User20Elmt.h)
set(src ${src} src/ElmtSystem/User20Elmt.cpp)

#############################################################
### For Materials system in AsFem                         ###
#############################################################
set(inc ${inc} include/MateSystem/MateType.h)
set(inc ${inc} include/MateSystem/MateNameDefine.h)
set(inc ${inc} include/MateSystem/MateBlock.h)
set(inc ${inc} include/MateSystem/MateSystem.h)
set(src ${src} src/MateSystem/MateSystem.cpp)
### for complex materials class
set(inc ${inc} include/MateSystem/Materials.h)
set(src ${src} src/MateSystem/Materials.cpp)
### for bulk(base) materials
set(inc ${inc} include/MateSystem/BulkMaterialBase.h)
set(inc ${inc} include/MateSystem/MechanicsMaterialBase.h)
set(inc ${inc} include/MateSystem/PhaseFieldFractureMaterialBase.h)
set(inc ${inc} include/MateSystem/PlasticMaterialBase.h)
set(inc ${inc} include/MateSystem/MultiphycisMechanicsMaterialBase.h)
### For constant poisson material
set(inc ${inc} include/MateSystem/ConstPoissonMaterial.h)
set(src ${src} src/MateSystem/ConstPoissonMaterial.cpp)
### For constant diffusion material
set(inc ${inc} include/MateSystem/ConstDiffusionMaterial.h)
set(src ${src} src/MateSystem/ConstDiffusionMaterial.cpp)
# for elastic material
set(inc ${inc} include/MateSystem/LinearElasticMaterial.h)
set(src ${src} src/MateSystem/LinearElasticMaterial.cpp)
### For incremental small strain material
set(inc ${inc} include/MateSystem/IncrementSmallStrainMaterial.h)
set(src ${src} src/MateSystem/IncrementSmallStrainMaterial.cpp)
### for SaintVenant Material
set(inc ${inc} include/MateSystem/SaintVenantMaterial.h)
set(src ${src} src/MateSystem/SaintVenantMaterial.cpp)
### For NeoHookean material
set(inc ${inc} include/MateSystem/NeoHookeanMaterial.h)
set(src ${src} src/MateSystem/NeoHookeanMaterial.cpp)
### For 1D plastic material with linear hardening
set(inc ${inc} include/MateSystem/Plastic1DMaterial.h)
set(src ${src} src/MateSystem/Plastic1DMaterial.cpp)
### For J2 Plasticity material
set(inc ${inc} include/MateSystem/J2PlasticityMaterial.h)
set(inc ${inc} src/MateSystem/J2PlasticityMaterial.cpp)
### For free energy based materials
set(inc ${inc} include/MateSystem/FreeEnergyMaterialBase.h)
#   for double well free energy materials
set(inc ${inc} include/MateSystem/DoubleWellFreeEnergyMaterial.h)
set(src ${src} src/MateSystem/DoubleWellFreeEnergyMaterial.cpp)
#   for ideal solution free energy
set(inc ${inc} include/MateSystem/IdealSolutionFreeEnergyMaterial.h)
set(src ${src} src/MateSystem/IdealSolutionFreeEnergyMaterial.cpp)
### For Miehe's linear elastic phase field fracture model
set(inc ${inc} include/MateSystem/MieheFractureMaterial.h)
set(src ${src} src/MateSystem/MieheFractureMaterial.cpp)
### For diffusion coupled phase field fracture model
set(inc ${inc} include/MateSystem/DiffusionFractureMaterial.h)
set(src ${src} src/MateSystem/DiffusionFractureMaterial.cpp)
### for stress decomposition phase field fracture material
set(inc ${inc} include/MateSystem/StressDecompositionMaterial.h)
set(src ${src} src/MateSystem/StressDecompositionMaterial.cpp)
### for neohookean phase-field fracture material
set(inc ${inc} include/MateSystem/NeoHookeanPFFractureMaterial.h)
set(src ${src} src/MateSystem/NeoHookeanPFFractureMaterial.cpp)
### for kobayashi material
set(inc ${inc} include/MateSystem/KobayashiMaterial.h)
set(src ${src} src/MateSystem/KobayashiMaterial.cpp)
### for diffusion coupled neohookean
set(inc ${inc} include/MateSystem/DiffNeoHookeanMaterial.h)
set(src ${src} src/MateSystem/DiffNeoHookeanMaterial.cpp)
### for linear elastic coupled cahn-hilliard material
set(inc ${inc} include/MateSystem/LinearElasticCHMaterial.h)
set(src ${src} src/MateSystem/LinearElasticCHMaterial.cpp)
### for wave material
set(inc ${inc} include/MateSystem/WaveMaterial.h)
set(src ${src} src/MateSystem/WaveMaterial.cpp)
### for thermal material
set(inc ${inc} include/MateSystem/ThermalMaterial.h)
set(src ${src} src/MateSystem/ThermalMaterial.cpp)
### For UMAT1
set(inc ${inc} include/MateSystem/User1Material.h)
set(src ${src} src/MateSystem/User1Material.cpp)
### For UMAT2
set(inc ${inc} include/MateSystem/User2Material.h)
set(src ${src} src/MateSystem/User2Material.cpp)
### For UMAT3
set(inc ${inc} include/MateSystem/User3Material.h)
set(src ${src} src/MateSystem/User3Material.cpp)
### For UMAT4
set(inc ${inc} include/MateSystem/User4Material.h)
set(src ${src} src/MateSystem/User4Material.cpp)
### For UMAT5
set(inc ${inc} include/MateSystem/User5Material.h)
set(src ${src} src/MateSystem/User5Material.cpp)
### For UMAT6
set(inc ${inc} include/MateSystem/User6Material.h)
set(src ${src} src/MateSystem/User6Material.cpp)
### For UMAT7
set(inc ${inc} include/MateSystem/User7Material.h)
set(src ${src} src/MateSystem/User7Material.cpp)
### For UMAT8
set(inc ${inc} include/MateSystem/User8Material.h)
set(src ${src} src/MateSystem/User8Material.cpp)
### For UMAT9
set(inc ${inc} include/MateSystem/User9Material.h)
set(src ${src} src/MateSystem/User9Material.cpp)
### For UMAT10
set(inc ${inc} include/MateSystem/User10Material.h)
set(src ${src} src/MateSystem/User10Material.cpp)
###
set(inc ${inc} include/MateSystem/BulkMateSystem.h)
set(src ${src} src/MateSystem/BulkMateSystem.cpp)
set(src ${src} src/MateSystem/InitBulkMateLibs.cpp)
set(src ${src} src/MateSystem/RunBulkMateLibs.cpp)
### for UMAT

#############################################################
### For Boundary condition system in AsFem                ###
#############################################################
set(inc ${inc} include/BCSystem/BCType.h)
set(inc ${inc} include/BCSystem/BCBlock.h)
### For integrated boundary conditions
set(inc ${inc} include/BCSystem/IntegrateBCBase.h)
### for neumann bc
set(inc ${inc} include/BCSystem/NeumannBC.h)
set(src ${src} src/BCSystem/NeumannBC.cpp)
### for flux bc
set(inc ${inc} include/BCSystem/FluxBC.h)
set(src ${src} src/BCSystem/FluxBC.cpp)
### for nodalbc base
set(inc ${inc} include/BCSystem/NodalBCBase.h)
### for nodal force bc
set(inc ${inc} include/BCSystem/NodalForceBC.h)
set(src ${src} src/BCSystem/NodalForceBC.cpp)
### for nodal flux bc
set(inc ${inc} include/BCSystem/NodalFluxBC.h)
set(src ${src} src/BCSystem/NodalFluxBC.cpp)
### for dirichlet bc base class
set(inc ${inc} include/BCSystem/DirichletBCBase.h)
### for dirichlet boundary condition
set(inc ${inc} include/BCSystem/DirichletBC.h)
set(src ${src} src/BCSystem/DirichletBC.cpp)
### for cyclic dirichlet bc
set(inc ${inc} include/BCSystem/CyclicDirichletBC.h)
set(src ${src} src/BCSystem/CyclicDirichletBC.cpp)
### for user1 dirichlet bc
set(inc ${inc} include/BCSystem/User1DirichletBC.h)
set(src ${src} src/BCSystem/User1DirichletBC.cpp)
### for user2 dirichlet bc
set(inc ${inc} include/BCSystem/User2DirichletBC.h)
set(src ${src} src/BCSystem/User2DirichletBC.cpp)
### for user3 dirichlet bc
set(inc ${inc} include/BCSystem/User3DirichletBC.h)
set(src ${src} src/BCSystem/User3DirichletBC.cpp)
### for user4 dirichlet bc
set(inc ${inc} include/BCSystem/User4DirichletBC.h)
set(src ${src} src/BCSystem/User4DirichletBC.cpp)
### for user5 dirichlet bc
set(inc ${inc} include/BCSystem/User5DirichletBC.h)
set(src ${src} src/BCSystem/User5DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/BCSystem.h)
###
set(src ${src} src/BCSystem/BCSystem.cpp)
set(src ${src} src/BCSystem/CheckAppliedBCNameIsValid.cpp)
set(src ${src} src/BCSystem/PrintBCSystemInfo.cpp)
set(src ${src} src/BCSystem/ApplyBC.cpp)
set(src ${src} src/BCSystem/ApplyDirichletBC.cpp)
set(src ${src} src/BCSystem/ApplyPresetBC.cpp)
set(src ${src} src/BCSystem/ApplyNodalDirichletBC.cpp)
set(src ${src} src/BCSystem/ApplyNodalNeumannBC.cpp)
set(src ${src} src/BCSystem/RunBCLibs.cpp)
### for user-defined-integrated-type bc
### for user1-bc
set(inc ${inc} include/BCSystem/User1BC.h)
set(src ${src} src/BCSystem/User1BC.cpp)
### for user2-bc
set(inc ${inc} include/BCSystem/User2BC.h)
set(src ${src} src/BCSystem/User2BC.cpp) 
### for user3-bc  
set(inc ${inc} include/BCSystem/User3BC.h)   
set(src ${src} src/BCSystem/User3BC.cpp)  
### for user4-bc  
set(inc ${inc} include/BCSystem/User4BC.h)
set(src ${src} src/BCSystem/User4BC.cpp)
### for user5-bc
set(inc ${inc} include/BCSystem/User5BC.h)
set(src ${src} src/BCSystem/User5BC.cpp)
### for user6-bc
set(inc ${inc} include/BCSystem/User6BC.h)
set(src ${src} src/BCSystem/User6BC.cpp)
### for user7-bc
set(inc ${inc} include/BCSystem/User7BC.h)
set(src ${src} src/BCSystem/User7BC.cpp)
### for user8-bc
set(inc ${inc} include/BCSystem/User8BC.h)
set(src ${src} src/BCSystem/User8BC.cpp)
### for user9-bc
set(inc ${inc} include/BCSystem/User9BC.h)
set(src ${src} src/BCSystem/User9BC.cpp)
### for user10-bc
set(inc ${inc} include/BCSystem/User10BC.h)
set(src ${src} src/BCSystem/User10BC.cpp)

#############################################################
### For Initial condition system in AsFem                 ###
#############################################################
set(inc ${inc} include/ICSystem/ICType.h)
set(inc ${inc} include/ICSystem/ICBlock.h)
set(inc ${inc} include/ICSystem/ICSystem.h)
set(src ${src} src/ICSystem/ICSystem.cpp)
set(src ${src} src/ICSystem/PrintICSystemInfo.cpp)
set(src ${src} src/ICSystem/ApplyIC.cpp)
set(src ${src} src/ICSystem/RunICLibs.cpp)
set(src ${src} src/ICSystem/ApplyConstantIC.cpp)
set(src ${src} src/ICSystem/ApplyRandomIC.cpp)
set(src ${src} src/ICSystem/ApplyRectangleIC.cpp)
set(src ${src} src/ICSystem/ApplyCircleIC.cpp)
set(src ${src} src/ICSystem/ApplySmoothCirlceIC.cpp)
set(src ${src} src/ICSystem/ApplyCubicIC.cpp)
set(src ${src} src/ICSystem/ApplySphereIC.cpp)
set(src ${src} src/ICSystem/User1IC.cpp)

#############################################################
### For FE space and the related functions in AsFem       ###
#############################################################
### for gauss-point
#########################
set(inc ${inc} include/FE/QPointType.h include/FE/QPointBase.h)
set(inc ${inc} include/FE/QPointGaussLegendre.h)
set(src ${src} src/FE/QPoints/QPointGaussLegendre.cpp)
set(src ${src} src/FE/QPoints/CreateGaussLegendrePoint.cpp)
set(inc ${inc} include/FE/QPointGaussLobatto.h)
set(src ${src} src/FE/QPoints/QPointGaussLobatto.cpp)
set(src ${src} src/FE/QPoints/CreateGaussLobattoPoint.cpp)
set(inc ${inc} include/FE/QPoint.h)
set(src ${src} src/FE/QPoints/QPoint.cpp)
set(src ${src} src/FE/QPoints/PrintQPointInfo.cpp)
#############################################################
### for FE
#########################
set(inc ${inc} include/FE/FE.h)
set(src ${src} src/FE/FE.cpp)
### for shape functions
set(inc ${inc} include/FE/LagrangeShapeFun.h)
### for 1D lagrange shape function
set(inc ${inc} include/FE/Lagrange1DShapeFun.h)
set(src ${src} src/FE/ShapeFuns/Lagrange1DShapeFun.cpp)
### for 2D lagrange shape function 
set(src ${src} include/FE/Lagrange2DShapeFun.h)
set(src ${src} src/FE/ShapeFuns/Lagrange2DShapeFun.cpp)
### for 3D lagrange shape function
set(inc ${inc} include/FE/Lagrange3DShapeFun.h)
set(src ${src} src/FE/ShapeFuns/Lagrange3DShapeFun.cpp)
###
set(src ${src} src/FE/ShapeFuns/LagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFuns/LagrangeShapeFunCalc.cpp)
set(inc ${inc} include/FE/ShapeFun.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun.cpp)
#############################################################
### For solution system in AsFem                          ###
#############################################################
set(inc ${inc} include/SolutionSystem/SolutionSystem.h)
set(src ${src} src/SolutionSystem/SolutionSystem.cpp)
set(src ${src} src/SolutionSystem/InitSolution.cpp)
set(src ${src} src/SolutionSystem/UpdateMaterials.cpp)

#############################################################
### For equation system in AsFem                          ###
#############################################################
set(inc ${inc} include/EquationSystem/EquationSystem.h)
set(src ${src} src/EquationSystem/EquationSystem.cpp)
set(src ${src} src/EquationSystem/CreateSparsityPattern.cpp)

#############################################################
### For nonlinear solver system in AsFem                  ###
#############################################################
set(inc ${inc} include/NonlinearSolver/NonlinearSolverType.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolverBlock.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolver.h)
set(src ${src} src/NonlinearSolver/NonlinearSolver.cpp)
set(src ${src} src/NonlinearSolver/Solve.cpp)
set(src ${src} src/NonlinearSolver/NewtowRaphsonSolve.cpp)

#############################################################
### For time stepping system in AsFem                     ###
#############################################################
set(inc ${inc} include/TimeStepping/TimeSteppingType.h)
set(inc ${inc} include/TimeStepping/TimeSteppingBlock.h)
set(inc ${inc} include/TimeStepping/TimeStepping.h)
set(src ${src} src/TimeStepping/TimeStepping.cpp)
set(src ${src} src/TimeStepping/Solve.cpp)

#############################################################
### For output system in AsFem                            ###
#############################################################
set(inc ${inc} include/OutputSystem/OutputType.h)
set(inc ${inc} include/OutputSystem/OutputBlock.h)
set(inc ${inc} include/OutputSystem/OutputSystem.h)
set(src ${src} src/OutputSystem/OutputSystem.cpp)
set(src ${src} src/OutputSystem/WriteResultToFile.cpp)
set(src ${src} src/OutputSystem/WriteResult2VTU.cpp)
set(src ${src} src/OutputSystem/WritePVDFile.cpp)

#############################################################
### For FE system in AsFem                                ###
#############################################################
set(inc ${inc} include/FESystem/FECalcType.h)
set(inc ${inc} include/FESystem/FESystem.h)
set(src ${src} src/FESystem/FESystem.cpp src/FESystem/InitBulkFESystem.cpp)
set(src ${src} src/FESystem/FormBulkFE.cpp)
set(src ${src} src/FESystem/FEAssemble.cpp)
set(src ${src} src/FESystem/FEProjection.cpp)

#############################################################
### For postprocess system in AsFem                       ###
#############################################################
set(inc ${inc} include/Postprocess/PostprocessType.h)
set(inc ${inc} include/Postprocess/PostprocessBlock.h)
set(inc ${inc} include/Postprocess/Postprocess.h)
set(src ${src} src/Postprocess/Postprocess.cpp)
### for settings check
set(src ${src} src/Postprocess/CheckWhetherPPSIsValid.cpp)
set(src ${src} src/Postprocess/RunPostprocess.cpp)
### for node pps
set(src ${src} src/Postprocess/NodeValuePostProcess.cpp)
### for elmt pps
set(src ${src} src/Postprocess/ElementValuePostProcess.cpp)
set(src ${src} src/Postprocess/AreaPostProcess.cpp)
set(src ${src} src/Postprocess/VolumePostProcess.cpp)
set(src ${src} src/Postprocess/ElementalIntegralPostProcess.cpp)
### for side-integral of projected variables
set(src ${src} src/Postprocess/SideIntegralPostProcess.cpp)
set(src ${src} src/Postprocess/ProjVariableSideIntegralPostProcess.cpp)
set(src ${src} src/Postprocess/Rank2MateSideIntegralPostProcess.cpp)

#############################################################
### For FEProblem in AsFem  (top level)                   ###
#############################################################
set(inc ${inc} include/FEProblem/FEControlInfo.h)
set(inc ${inc} include/FEProblem/FEJobType.h)
set(inc ${inc} include/FEProblem/FEJobBlock.h)
set(inc ${inc} include/FEProblem/FEProblem.h)
set(src ${src} src/FEProblem/FEProblem.cpp)
set(src ${src} src/FEProblem/PreRun.cpp)
set(src ${src} src/FEProblem/Run.cpp)
set(src ${src} src/FEProblem/RunStaticAnalysis.cpp)
set(src ${src} src/FEProblem/RunTransientAnalysis.cpp)

##################################################
add_executable(asfem ${inc} ${src})


##################################################
### Following lines are used by vim            ###
### you can delete all of them                 ###
##################################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PETSC_DIR}/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_DIR}/include")

