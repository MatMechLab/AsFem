cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 17)

###############################################
### Set your PETSc path here                ###
###############################################
set(PETSC_DIR "/home/by/Programs/petsc/3.13.0-opt")
set(MPI_DIR "/home/by/Programs/mpich/3.3.2")
set(Eigen "/home/by/Programs/Eigen-dev")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


# For PETSc
include_directories("${PETSC_DIR}/include")
link_libraries("${PETSC_DIR}/lib/libpetsc.so")
# For MPI
include_directories("${MPI_DIR}/include")
link_libraries("${MPI_DIR}/lib/libmpi.so")
# For Eigen
include_directories("${Eigen}")

###############################################
### For debug mode                          ###
###############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O2 -fopenmp")
### for gcc 9, you must add this!!!
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=deprecated-copy")


###############################################
### For release mode                        ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG")
# if you have mumps installed, please use the following line
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG -DHASMUMPS")


###############################################
### For release mode in Windows (VS2019)    ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W1 /arch:AVX")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /openmp")


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
### For Utils for string                                  ###
#############################################################
set(inc ${inc} include/Utils/StringUtils.h)
set(src ${src} src/Utils/StringUtils/StringCaseConvert.cpp)
set(src ${src} src/Utils/StringUtils/RemoveString.cpp)
set(src ${src} src/Utils/StringUtils/GoToLine.cpp)
set(src ${src} src/Utils/StringUtils/IsBracketMatch.cpp)
set(src ${src} src/Utils/StringUtils/IsCommentLine.cpp)
set(src ${src} src/Utils/StringUtils/IsUniqueStrVec.cpp)
set(src ${src} src/Utils/StringUtils/IsValidExpression.cpp)
set(src ${src} src/Utils/StringUtils/SplitStr.cpp)
set(src ${src} src/Utils/StringUtils/SplitStrNum.cpp)
set(src ${src} src/Utils/StringUtils/SplitStrNumAfter.cpp)

#############################################################
### For Math Utils function                               ###
#############################################################
set(inc ${inc} include/Utils/MathUtils.h)
set(inc ${inc} include/Utils/Vector3d.h)
set(inc ${inc} include/Utils/MatrixXd.h)
set(inc ${inc} include/Utils/VectorXd.h)
### For Vector3d
set(src ${src} src/Utils/Vector3d.cpp)
### For rank-2 tensor
set(inc ${inc} include/Utils/RankTwoTensor.h)
set(src ${src} src/Utils/Tensors/RankTwoTensor.cpp)
### For rank-4 tensor
set(inc ${inc} include/Utils/RankFourTensor.h)
set(src ${src} src/Utils/Tensors/RankFourTensor.cpp)

#############################################################
### For message printer                                   ###
#############################################################
set(inc ${inc} include/MessagePrinter/MessagePrinter.h)
set(src ${src} src/MessagePrinter/Program/Msg_AsFem_Exit.cpp)
### For input file
set(src ${src} src/MessagePrinter/Input/Msg_Input_InvalidArgs.cpp)
set(src ${src} src/MessagePrinter/Input/Msg_Input_InvalidInputFileName.cpp)
set(src ${src} src/MessagePrinter/Input/Msg_Input_InvalidArgs.cpp)
set(src ${src} src/MessagePrinter/Input/Msg_Input_LineError.cpp)
set(src ${src} src/MessagePrinter/Input/Msg_InputFileName_Invalid.cpp)
set(src ${src} src/MessagePrinter/Input/Msg_Input_BlockBracketNotComplete.cpp)


#############################################################
### For InputSystem                                       ###
#############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp)
set(src ${src} src/InputSystem/ReadInputFile.cpp)
set(src ${src} src/InputSystem/ReadMeshBlock.cpp)
set(src ${src} src/InputSystem/ReadDofsBlock.cpp)
set(src ${src} src/InputSystem/ReadElmtBlock.cpp)
set(src ${src} src/InputSystem/ReadMateBlock.cpp)
set(src ${src} src/InputSystem/ReadBCBlock.cpp)
set(src ${src} src/InputSystem/ReadICBlock.cpp)
set(src ${src} src/InputSystem/ReadQPointBlock.cpp)
set(src ${src} src/InputSystem/ReadJobBlock.cpp)
set(src ${src} src/InputSystem/ReadProjectionBlock.cpp)
set(src ${src} src/InputSystem/ReadNonlinearSolverBlock.cpp)
set(src ${src} src/InputSystem/ReadTimeSteppingBlock.cpp)
set(src ${src} src/InputSystem/ReadOutputBlock.cpp)

#############################################################
### For Mesh class                                        ###
#############################################################
set(inc ${inc} include/Mesh/Mesh.h include/Mesh/MeshType.h)
set(inc ${inc} include/Mesh/Nodes.h)
set(src ${src} src/Mesh/Mesh.cpp src/Mesh/CreateMesh.cpp)
set(src ${src} src/Mesh/Create1DMesh.cpp src/Mesh/Create2DMesh.cpp src/Mesh/Create3DMesh.cpp)
set(src ${src} src/Mesh/PrintMeshInfo.cpp)
set(src ${src} src/Mesh/SaveMesh.cpp)
set(src ${src} src/Mesh/SetElmtInfo.cpp)
set(src ${src} src/Mesh/ReadMeshFromGmsh.cpp)
set(src ${src} src/Mesh/ReadMeshFromAbaqus.cpp)
set(src ${src} src/Mesh/GetGmshInfo.cpp)
set(src ${src} src/Mesh/GetAbaqusInfo.cpp)


#############################################################
### For DofHandler class                                  ###
#############################################################
set(inc ${inc} include/DofHandler/DofHandler.h)
set(src ${src} src/DofHandler/DofHandler.cpp)
set(src ${src} src/DofHandler/AddNameToDofNameList.cpp)
set(src ${src} src/DofHandler/CreateDofMap.cpp)
set(src ${src} src/DofHandler/PrintDofInfo.cpp)


#############################################################
### For Boundary system                                   ###
#############################################################
set(inc ${inc} include/BCs/BCBlock.h include/BCs/BCType.h)
set(inc ${inc} include/BCs/BCSystem.h)
set(src ${src} src/BCs/BCSystem.cpp src/BCs/AddBCBlock.cpp)
set(src ${src} src/BCs/InitBCSystem.cpp)
set(src ${src} src/BCs/PrintBCSystemInfo.cpp)
set(src ${src} src/BCs/ApplyBC.cpp)
set(src ${src} src/BCs/ApplyDirichletBC.cpp)
set(src ${src} src/BCs/ApplyPresetBC.cpp)
set(src ${src} src/BCs/ApplyNeumannBC.cpp)
set(src ${src} src/BCs/ApplyPressureBC.cpp)
set(src ${src} src/BCs/ApplyUser1BC.cpp)



#############################################################
### For Initial condition system                          ###
#############################################################
set(inc ${inc} include/ICs/ICBlock.h include/ICs/ICType.h)
set(inc ${inc} include/ICs/ICSystem.h)
set(src ${src} src/ICs/ICSystem.cpp src/ICs/AddICBlock.cpp)
set(src ${src} src/ICs/PrintICSystemInfo.cpp)
set(src ${src} src/ICs/ApplyIC.cpp)
set(src ${src} src/ICs/ApplyConstIC.cpp)
set(src ${src} src/ICs/ApplyRandomIC.cpp)
set(src ${src} src/ICs/ApplyCircleIC.cpp)

#############################################################
### For Element  System                                   ###
#############################################################
set(inc ${inc} include/ElmtSystem/ElmtBlock.h include/ElmtSystem/ElmtType.h)
set(inc ${inc} include/ElmtSystem/ElmtSystem.h src/ElmtSystem/ElmtSystem.cpp)
set(src ${src} src/ElmtSystem/PrintElmtSystemInfo.cpp)
set(src ${src} src/ElmtSystem/AddElmtBlock.cpp)
set(src ${src} src/ElmtSystem/RunElmtLib.cpp)
set(src ${src} src/ElmtSystem/Poisson.cpp)
set(src ${src} src/ElmtSystem/Mechanics.cpp)
set(src ${src} src/ElmtSystem/ThermalMechanics.cpp)
set(src ${src} src/ElmtSystem/Diffusion.cpp)
### for CahnHilliard elmts
set(src ${src} src/ElmtSystem/CahnHilliard.cpp)
set(src ${src} src/ElmtSystem/MechanicalCahnHilliard.cpp)
set(src ${src} src/ElmtSystem/TensorCahnHilliard.cpp)
set(src ${src} src/ElmtSystem/MieheFracture.cpp)
set(src ${src} src/ElmtSystem/CohesivePFFracture.cpp)
set(src ${src} src/ElmtSystem/BordenFracture.cpp)
### for thermal conduct problem
set(src ${src} src/ElmtSystem/Thermal.cpp)
### for dendrite model
set(src ${src} src/ElmtSystem/DendriteModel.cpp)
### for User-Element (UEL)
set(src ${src} src/ElmtSystem/UserElmt1.cpp)



#############################################################
### For Material System                                   ###
#############################################################
set(inc ${inc} include/MateSystem/MateBlock.h include/MateSystem/MateType.h)
set(inc ${inc} include/MateSystem/MateSystem.h)
set(src ${src} src/MateSystem/MateSystem.cpp)
set(src ${src} src/MateSystem/AddMateBlock.cpp)
set(src ${src} src/MateSystem/RunMateLib.cpp)
set(src ${src} src/MateSystem/InitMateSystem.cpp)
set(src ${src} src/MateSystem/PrintMateSystemInfo.cpp)
### for const poisson material
set(src ${src} src/MateSystem/ConstPoissonMaterial.cpp)
### for const poisson material
set(src ${src} src/MateSystem/NonlinearPoissonMaterial.cpp)
### for const diffusion material
set(src ${src} src/MateSystem/ConstDiffusionMaterial.cpp)
### for cahn-hilliard material
set(src ${src} src/MateSystem/CahnHilliardMaterial.cpp)
set(src ${src} src/MateSystem/TensorCahnHilliardMaterial.cpp)
### for mechanically coupled CahnHilliard material
set(src ${src} src/MateSystem/LinearElasticCahnHilliardMaterial.cpp)
### for linear-thermal-elastic-material
set(src ${src} src/MateSystem/LinearThermalElasticMaterial.cpp)


### for linear elastic material
set(src ${src} src/MateSystem/LinearElasticMaterial.cpp)
### for saint venant material
set(src ${src} src/MateSystem/SaintVenantMaterial.cpp)
### for compressive neohookean material
set(src ${src} src/MateSystem/NeoHookeanMaterial.cpp)

### for phase field fracture material 
set(src ${src} src/MateSystem/MieheLinearElasticMaterial.cpp)
### for phase field fracture with neohookean material
set(src ${src} src/MateSystem/MieheNeoHookeanMaterial.cpp)
### for anisotropic phase field fracture model
set(src ${src} src/MateSystem/AnisoPFFractureMaterial.cpp)
### for borden's phase field fracture material 
set(src ${src} src/MateSystem/BordenLinearElasticMaterial.cpp)
set(src ${src} src/MateSystem/CohesivePFFractureMaterial.cpp)
### for user's materials
set(src ${src} src/MateSystem/UserMaterial1.cpp)
set(src ${src} src/MateSystem/UserMaterial2.cpp)
### for dendrite material (in 2d)
set(src ${src} src/MateSystem/DendriteMaterial.cpp)
### for current induced thermal analysis
set(src ${src} src/MateSystem/CurrentThermalMaterial.cpp)


#############################################################
### For Solution System                                   ###
#############################################################
set(inc ${inc} include/Solution/Solution.h)
set(src ${src} src/Solution/Solution.cpp src/Solution/InitSolution.cpp)

#############################################################
### For Equation System                                   ###
#############################################################
set(inc ${inc} include/EquationSystem/EquationSystem.h)
set(src ${src} src/EquationSystem/EquationSystem.cpp)
set(src ${src} src/EquationSystem/InitEquationSystem.cpp)
set(src ${src} src/EquationSystem/CreateSparsityPattern.cpp)



#############################################################
### For FE space(shape fun and gauss integration)         ###
#############################################################
set(inc ${inc} include/FE/FE.h)
set(src ${src} src/FE/FE.cpp src/FE/InitFE.cpp)
set(src ${src} src/FE/PrintQPointInfo.cpp)
##################
### For ShapeFun
##################
set(inc ${inc} include/FE/ShapeFun.h)
set(src ${src} src/FE/ShapeFun/ShapeFun.cpp)
set(src ${src} src/FE/ShapeFun/ShapeFunCalc.cpp)
set(src ${src} src/FE/ShapeFun/Compute1DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFun/Compute2DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFun/Compute3DLagrangeShapeFun.cpp)

##################
### For QPoint
##################
set(inc ${inc} include/FE/QPoint.h)
set(src ${src} src/FE/QPoint/QPoint.cpp)
set(src ${src} src/FE/QPoint/CreateQPoints.cpp)
set(src ${src} src/FE/QPoint/Create1DGaussPoint.cpp)
set(src ${src} src/FE/QPoint/Create2DGaussPoint.cpp)
set(src ${src} src/FE/QPoint/Create3DGaussPoint.cpp)


#############################################################
### For FE System (form K and Residual...)                ###
#############################################################
set(inc ${inc} include/FESystem/FESystem.h)
set(src ${src} src/FESystem/FESystem.cpp)
set(src ${src} src/FESystem/InitFESystem.cpp)
set(src ${src} src/FESystem/FormFE.cpp)
set(src ${src} src/FESystem/AssembleLocalToGlobal.cpp)
set(src ${src} src/FESystem/AssembleLocalHistToGlobal.cpp)
set(src ${src} src/FESystem/Projection.cpp)

#############################################################
### For Nonlinear solver settings                         ###
#############################################################
set(inc ${inc} include/NonlinearSolver/NonlinearSolverBlock.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolverType.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolver.h)
set(src ${src} src/NonlinearSolver/NonlinearSolver.cpp)
set(src ${src} src/NonlinearSolver/Init.cpp)
set(src ${src} src/NonlinearSolver/SSolve.cpp)
set(src ${src} src/NonlinearSolver/Solve.cpp)
set(src ${src} src/NonlinearSolver/LinearSolve.cpp)
set(src ${src} src/NonlinearSolver/NewtonRaphson.cpp)
set(src ${src} src/NonlinearSolver/CheckConvergence.cpp)
set(src ${src} src/NonlinearSolver/PrintIterationInfo.cpp)


#############################################################
### For Output System                                     ###
#############################################################
set(inc ${inc} include/OutputSystem/OutputBlock.h include/OutputSystem/OutputFileType.h)
set(inc ${inc} include/OutputSystem/OutputSystem.h)
set(src ${src} src/OutputSystem/OutputSystem.cpp)
set(src ${src} src/OutputSystem/InitOutputStream.cpp)
set(src ${src} src/OutputSystem/WriteResultToVTU.cpp)


#############################################################
### For Time stepping                                     ###
#############################################################
set(inc ${inc} include/TimeStepping/TimeSteppingBlock.h)
set(inc ${inc} include/TimeStepping/TimeSteppingType.h)
set(inc ${inc} include/TimeStepping/TimeStepping.h)
set(src ${src} src/TimeStepping/TimeStepping.cpp)
set(src ${src} src/TimeStepping/Init.cpp)
set(src ${src} src/TimeStepping/Stepping.cpp)
set(src ${src} src/TimeStepping/SteppingNew.cpp)
### for backward euler
set(src ${src} src/TimeStepping/BackwardEuler.cpp)
### for cranck-nicolson method
set(src ${src} src/TimeStepping/CrankNicolson.cpp)

#############################################################
### For postprocess system                                ###
#############################################################
set(inc ${inc} include/PPSSystem/PPSSystem.h)


#############################################################
### For FE Problem(the top level of AsFem)                ###
#############################################################
### For Job Block
set(inc ${inc} include/FEProblem/JobBlock.h include/FEProblem/JobType.h)
set(inc ${inc} include/FEProblem/FeCtrlInfo.h)
set(inc ${inc} include/FEProblem/FEProblem.h)
set(src ${src} src/FEProblem/FEProblem.cpp)
set(src ${src} src/FEProblem/RunFEProblem.cpp)
set(src ${src} src/FEProblem/RunInputSystem.cpp)
set(src ${src} src/FEProblem/InitFEProblem.cpp)
set(src ${src} src/FEProblem/RunAnalysis.cpp)
set(src ${src} src/FEProblem/RunStaticAnalysis.cpp)
set(src ${src} src/FEProblem/RunTransientAnalysis.cpp)


##################################################
add_executable(asfem ${inc} ${src})

