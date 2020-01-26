cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 14)

###############################################
### Set your Eigen path here                ###
###############################################
set(Eigen "/home/by/Programs/Eigen-dev")

#set(SuiteSparse "/home/by/Programs/SuiteSparse")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${Eigen})

#include_directories("${SuiteSparse}/include")
#link_libraries("${SuiteSparse}/lib/libumfpack.so")

###############################################
### For debug mode                          ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -fopenmp -O2 -DUSE_SUITESPARSE")

###############################################
### For release mode                        ###
###############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG")

###############################################
### For windows release                     ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W1 /favor:INTEL64 /arch:AVX")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GL /Qpar /openmp")
###################################################

###############################################
### For intel compiler                      ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -qopenmp -mkl:parallel -inline-forceinline")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -DUSE_PARDISO")

# the following line takes too long time for compiling
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xAVX -axCORE-AVX2 -DUSE_PARDISO")

#############################################################
#############################################################
### For beginners, please don't edit the following line!  ###
### Do not edit the following lines, please!!!            ###
### Do not edit the following lines, please!!!            ###
### Do not edit the following lines, please!!!            ###
#############################################################
#############################################################
# For Welcome header file and main.cpp
set(inc include/Welcome.h)
set(src src/main.cpp)


#############################################################
### For message print system                              ###
#############################################################
set(inc ${inc} include/MessagePrint/MessagePrint.h)
set(src ${src} src/MessagePrint/Program/Msg_AsFem_Exit.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_InvalidArgs.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFileName_Invalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LineError.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoMeshBlock.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_InvalidMeshType.cpp)
### for [dofs] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_DofNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_DofNameDuplicate.cpp)
### for [projection] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ProjectionNameNotFound.cpp)
### for [bcs] block
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoDofsBlock.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BlockBracketNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCSubBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCUserNumNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCUserNumNotValid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NoDofNameFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCBoundaryNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCValueNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCBlockInfoNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NoBCBlockWarning.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_BCBlockUnsupportedType.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_BCBlockBlockNameInvalid.cpp)
### for [qpoint] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_QpBlockTypeNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_QpBlockOrderNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_QpBlockOrderInvalid.cpp)
### for [linearsolver] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LinearSolverBlock_InvalidSolverOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LinearSolverBlock_MaxIterNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LinearSolverBlock_MaxIterInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LinearSolverBlock_TolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_LinearSolverBlock_TolInvalid.cpp)
### For [nonlinearsolver] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_InvalidSolverOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_MaxIterNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_MaxIterInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_RrelTolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_RrelTolInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_RabsTolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_RabsTolInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_ErelTolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_ErelTolInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_EabsTolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_EabsTolInvalid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_STolNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_NonLinearSolverBlock_STolInvalid.cpp)



### for [elmts] block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtUserNumNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtUserNumNotValid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtBlockUnsupportedType.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtBlockMateNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtBlockBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtBlockInfoNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtSubBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ElmtBlockDuplicateDofsName.cpp)

set(src ${src} src/MessagePrint/InputSystem/Msg_ElmtBlockBlockNameInvalid.cpp)



set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoElmtBlockFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICSubBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICUserNumNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICUserNumNotValid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICValueNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICBlockInfoNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICBlockUnsupportedType.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoICBlockFoundWarning.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_ICBlockNoValueWarning.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_ICBlockBlockNameInvalid.cpp)

### for input check
set(src ${src} src/MessagePrint/Program/Msg_InputCheck_DofInvalidInElmtBlock.cpp)
set(src ${src} src/MessagePrint/Program/Msg_InputCheck_DofInvalidInBCBlock.cpp)
set(src ${src} src/MessagePrint/Program/Msg_InputCheck_DofInvalidInICBlock.cpp)
set(src ${src} src/MessagePrint/Program/Msg_InputCheck_ElmtBlockDuplicateDomain.cpp)
set(src ${src} src/MessagePrint/Program/Msg_InputCheck_DofIsNotAssignedToElmt.cpp)


set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateSubBlockNameNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateUserNumNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateUserNumNotValid.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateParamsNotFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateBlockInfoNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateBlockUnsupportType.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_MateBlockNoParamsWarning.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoMateBlockFoundWarning.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputCheck_MateNotFoundForElmtBlock.cpp)

### For job block
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidJobType.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidDebugOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidSolverOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidNonLinearSolverOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidTimeSteppingOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_InputFile_NoJobBlockFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_NoTimeFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidTimeValue.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_NoDtFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidDtValue.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBblock_TransientInfoNotComplete.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_NoIntervalFound.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidIntervalValue.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidProjectionOption.cpp)
set(src ${src} src/MessagePrint/InputSystem/Msg_Input_JobBlock_InvalidAdaptiveOption.cpp)


### For gauss point generation
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid1DLintNum.cpp)
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid1DLobbotLintNum.cpp)

set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid2DGaussOrder.cpp)
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid2DMeshType.cpp)
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid2DLobattoMeshType.cpp)

set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid3DGaussOrder.cpp)
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid3DMeshType.cpp)
set(src ${src} src/MessagePrint/GaussPoint/Msg_GaussPoint_Invalid3DLobattoMeshType.cpp)

set(src ${src} src/MessagePrint/ShapeFun/Msg_ShapeFun_InformationNotComlete.cpp)
set(src ${src} src/MessagePrint/ShapeFun/Msg_ShapeFun_UnsupportedMeshType.cpp)

set(src ${src} src/MessagePrint/ShapeFun/Msg_ShapeFun_Lagrange1DElmtSingular.cpp)
set(src ${src} src/MessagePrint/ShapeFun/Msg_ShapeFun_Lagrange2DElmtSingular.cpp)
set(src ${src} src/MessagePrint/ShapeFun/Msg_ShapeFun_Lagrange3DElmtSingular.cpp)

#### for linear solver
set(src ${src} src/MessagePrint/LinearSolver/Msg_LinearSolver_AnalyzeFailed.cpp)
set(src ${src} src/MessagePrint/LinearSolver/Msg_LinearSolver_FactorizeFailed.cpp)
set(src ${src} src/MessagePrint/LinearSolver/Msg_LinearSolver_SolveFailed.cpp)

#### for nonlinear solver
set(src ${src} src/MessagePrint/NonLinearSolver/Msg_NonLinearSolver_UnsupportSolver.cpp)


#### for equation system
set(src ${src} src/MessagePrint/EquationSystem/Msg_EquationSystem_NotInitialized.cpp)

#### for Material System
set(src ${src} src/MessagePrint/MaterialSystem/Msg_MaterialSystem_UnsupportUmat.cpp)

#############################################################
### For string utils                                      ###
#############################################################
set(inc ${inc} include/Utils/StringUtils.h)
set(src ${src} src/Utils/StringUtils/RemoveStrSpace.cpp)
set(src ${src} src/Utils/StringUtils/RemoveSymbolFromStr.cpp)
set(src ${src} src/Utils/StringUtils/SplitStrNum.cpp)
set(src ${src} src/Utils/StringUtils/SplitStrNumAfter.cpp)
set(src ${src} src/Utils/StringUtils/SplitStr.cpp)
set(src ${src} src/Utils/StringUtils/IsCommentLine.cpp)
set(src ${src} src/Utils/StringUtils/IsBracketMatch.cpp)
set(src ${src} src/Utils/StringUtils/GoToLine.cpp)
set(src ${src} src/Utils/StringUtils/IsUniqueStrVec.cpp)
set(src ${src} src/Utils/StringUtils/IsValidExpression.cpp)

#############################################################
### For input system                                      ###
#############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp src/InputSystem/InitInputSystem.cpp)
set(src ${src} src/InputSystem/ReadInputFile.cpp)
set(src ${src} src/InputSystem/ReadMeshBlock.cpp)
set(src ${src} src/InputSystem/ReadDofsBlock.cpp)
set(src ${src} src/InputSystem/ReadBCBlock.cpp)
set(src ${src} src/InputSystem/ReadElmtBlock.cpp)
set(src ${src} src/InputSystem/ReadICBlock.cpp)
set(src ${src} src/InputSystem/ReadMateBlock.cpp)
set(src ${src} src/InputSystem/ReadJobBlock.cpp)
set(src ${src} src/InputSystem/ReadQpBlock.cpp)
set(src ${src} src/InputSystem/ReadLinearSolverBlock.cpp)
set(src ${src} src/InputSystem/ReadNonLinearSolverBlock.cpp)
set(src ${src} src/InputSystem/ReadProjectionBlock.cpp)



#############################################################
### For mesh system                                       ###
#############################################################
set(inc ${inc} include/Mesh/Mesh.h include/Mesh/MeshTypeDefine.h)
set(inc ${inc} include/Mesh/Point.h include/Mesh/Nodes.h)
set(src ${src} src/Mesh/Mesh.cpp src/Mesh/CreateMesh.cpp src/Mesh/PrintMeshInfo.cpp)
set(src ${src} src/Mesh/Create1DMesh.cpp src/Mesh/Create2DMesh.cpp src/Mesh/Create3DMesh.cpp)
set(src ${src} src/Mesh/SetMeshElmtInfo.cpp)
set(src ${src} src/Mesh/ReadMeshFromGmsh.cpp)
set(src ${src} src/Mesh/SaveMesh.cpp)
set(src ${src} src/Mesh/GetGmshInfo.cpp)

#############################################################
### For dof handler system                                ###
#############################################################
set(inc ${inc} include/DofHandler/DofHandler.h)
set(src ${src} src/DofHandler/DofHandler.cpp src/DofHandler/DofHandlerSettings.cpp)
set(src ${src} src/DofHandler/CreateDofMap.cpp src/DofHandler/Init.cpp)
set(src ${src} src/DofHandler/DofHandlerPrinter.cpp)
set(src ${src} src/DofHandler/ModifyDofActiveMapViaBC.cpp)


#############################################################
### For boundary condition system                         ###
#############################################################
set(inc ${inc} include/BCs/BCSystem.h include/BCs/BCBlock.h)
set(src ${src} src/BCs/BCSystem.cpp src/BCs/AddBCBlock.cpp)
set(src ${src} src/BCs/Init.cpp)
set(src ${src} src/BCs/PrintBCBlockInfo.cpp)
set(src ${src} src/BCs/ApplyBC.cpp)
set(src ${src} src/BCs/ApplyPresetBC.cpp)
set(src ${src} src/BCs/ApplyDirichletBC.cpp)
set(src ${src} src/BCs/ApplyNeumannBC.cpp)
set(src ${src} src/BCs/ApplyPreSetDirichletBC.cpp)


#############################################################
### For initial condition system                          ###
#############################################################
set(inc ${inc} include/ICs/ICBlock.h include/ICs/ICType.h)
set(inc ${inc} include/ICs/ICSystem.h)
set(src ${src} src/ICs/ICSystem.cpp src/ICs/AddICBlock.cpp)
set(src ${src} src/ICs/PrintICBlockInfo.cpp)
set(src ${src} src/ICs/ApplyIC.cpp)
set(src ${src} src/ICs/ApplyConstIC.cpp)
set(src ${src} src/ICs/ApplyRandIC.cpp)
set(src ${src} src/ICs/ApplyCircleIC.cpp)


#############################################################
### For element system                                    ###
#############################################################
set(inc ${inc} include/ElmtSystem/ElmtBlock.h include/ElmtSystem/ElmtSystem.h)
SET(inc ${inc} include/ElmtSystem/ElmtType.h)
set(src ${src} src/ElmtSystem/ElmtSystem.cpp)
set(src ${src} src/ElmtSystem/AddElmtBlock.cpp)
set(src ${src} src/ElmtSystem/PrintElmtBlockInfo.cpp)
set(src ${src} src/ElmtSystem/RunElmtLib.cpp)
### for poisson equation
set(src ${src} src/ElmtSystem/Poisson.cpp)
### for diffusion equation
set(src ${src} src/ElmtSystem/Diffusion.cpp)
### for cahn hilliard equation
set(src ${src} src/ElmtSystem/CahnHilliard.cpp)
### for mechanics
set(src ${src} src/ElmtSystem/Mechanics.cpp)
### for wave equation
set(src ${src} src/ElmtSystem/Wave.cpp)
### for thermal mechanics
set(src ${src} src/ElmtSystem/ThermalMechanics.cpp)
### for linear elastic phase field fracture model
set(src ${src} src/ElmtSystem/LinearElasticPhaseFieldFracture.cpp)
### for user-defined element (UEL)
set(src ${src} src/ElmtSystem/UserElmt1.cpp)


#############################################################
### For material system                                   ###
#############################################################
set(inc ${inc} include/MaterialSystem/MaterialSystem.h include/MaterialSystem/MaterialBlock.h)
set(inc ${inc} include/MaterialSystem/MaterialType.h)
set(src ${src} src/MaterialSystem/MaterialSystem.cpp)
set(src ${src} src/MaterialSystem/AddMaterialBlock.cpp)
set(src ${src} src/MaterialSystem/PrintMateBlockInfo.cpp)
set(src ${src} src/MaterialSystem/RunMateLib.cpp)
set(src ${src} src/MaterialSystem/ConstPoissonMate.cpp)
set(src ${src} src/MaterialSystem/NonlinearPoissonMate.cpp)
set(src ${src} src/MaterialSystem/ConstDiffusionMate.cpp)
set(src ${src} src/MaterialSystem/NonlinearDiffusionMate.cpp)
### For CahnHilliard material
set(src ${src} src/MaterialSystem/CahnHilliardMate.cpp)
### For mechanics material
set(src ${src} src/MaterialSystem/LinearElasticMate.cpp)
### For SaintVenantMate
set(src ${src} src/MaterialSystem/SaintVenantMate.cpp)
### For finite strain material
set(src ${src} src/MaterialSystem/FiniteStrainMate.cpp)
### For NeoHookeanMate material
set(src ${src} src/MaterialSystem/NeoHookeanMate.cpp)
### For thermal elastic material
set(src ${src} src/MaterialSystem/ThermalElasticMate.cpp)

### For linear elastic phase field fracture model (From Miehe's model)
set(src ${src} src/MaterialSystem/LinearElasticPhaseFieldFracMate.cpp)
set(src ${src} src/MaterialSystem/ModifyElasticFracMate.cpp)


### For general wave equation's material
set(src ${src} src/MaterialSystem/GenralWaveMate.cpp)
### UMAT
set(src ${src} src/MaterialSystem/UserMate1.cpp)

#############################################################
### For finite element system                             ###
#############################################################
set(inc ${inc} include/FE/FE.h)
set(src ${src} src/FE/FE.cpp src/FE/InitFE.cpp)


### for qp points
set(inc ${inc} include/FE/QPoint.h include/FE/QPBlockInfo.h)
set(src ${src} src/FE/GaussPoint/QPoint.cpp)
set(src ${src} src/FE/GaussPoint/CreateQPoints.cpp)
set(src ${src} src/FE/GaussPoint/Create1DGaussPoint.cpp)
set(src ${src} src/FE/GaussPoint/Create2DGaussPoint.cpp)
set(src ${src} src/FE/GaussPoint/Create3DGaussPoint.cpp)

### for shape funs
set(inc ${inc} include/FE/ShapeFun.h)
set(src ${src} src/FE/ShapeFun/ShapeFun.cpp src/FE/ShapeFun/ShapeFunCalc.cpp)
set(src ${src} src/FE/ShapeFun/Compute1DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFun/Compute2DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFun/Compute3DLagrangeShapeFun.cpp)

#############################################################
### For equation system                                   ###
#############################################################
set(inc ${inc} include/EquationSystem/EquationSystem.h)
set(src ${src} src/EquationSystem/EquationSystem.cpp)
set(src ${src} src/EquationSystem/InitEquationSystem.cpp)
set(src ${src} src/EquationSystem/CreateSparsityPatterns.cpp)


#############################################################
### For solution system                                   ###
#############################################################
set(inc ${inc} include/SolutionSystem/SolutionSystem.h)
set(src ${src} src/SolutionSystem/SolutionSystem.cpp)
set(src ${src} src/SolutionSystem/SolutionSystemInit.cpp)
set(src ${src} src/SolutionSystem/UpdateSolution.cpp)


#############################################################
### For linear solver system                              ###
#############################################################
set(inc ${inc} include/LinearSolver/LinearSolver.h include/LinearSolver/LinearSolverType.h)
set(inc ${inc} include/LinearSolver/LinearSolverBlockInfo.h)
set(src ${src} src/LinearSolver/LinearSolver.cpp)
set(src ${src} src/LinearSolver/LinearSolverGetInfo.cpp)
set(src ${src} src/LinearSolver/Solve.cpp)
set(src ${src} src/LinearSolver/InitLinearSolver.cpp)
set(src ${src} src/LinearSolver/PrintSolverInfo.cpp)


#############################################################
### For nonlinear solver system                           ###
#############################################################
set(inc ${inc} include/NonLinearSolver/NonLinearSolver.h)
set(inc ${inc} include/NonLinearSolver/NonLinearSolverType.h)
set(inc ${inc} include/NonLinearSolver/NonLinearSolverBlockInfo.h)
set(src ${src} src/NonLinearSolver/NonLinearSolver.cpp)
set(src ${src} src/NonLinearSolver/PrintNonLinearSolverInfo.cpp)
set(src ${src} src/NonLinearSolver/InitNonlinearSolver.cpp)
set(src ${src} src/NonLinearSolver/NonLinearSolve.cpp)
set(src ${src} src/NonLinearSolver/ConvergenceCheck.cpp)
set(src ${src} src/NonLinearSolver/Linear.cpp)
set(src ${src} src/NonLinearSolver/NewtonRaphson.cpp)
set(src ${src} src/NonLinearSolver/NewtonRaphsonModify.cpp)
set(src ${src} src/NonLinearSolver/NewtonRaphsonWithLineSearch.cpp)
set(src ${src} src/NonLinearSolver/NewtonRaphsonWithBisectionLineSearch.cpp)
set(src ${src} src/NonLinearSolver/NewtonRaphsonWithRegulaFalsiLineSearch.cpp)
set(src ${src} src/NonLinearSolver/PrintIterationInfo.cpp)


#############################################################
### For time stepping method system                       ###
#############################################################
set(inc ${inc} include/TimeStepping/TimeStepping.h include/TimeStepping/TimeSteppingType.h)
set(src ${src} src/TimeStepping/TimeStepping.cpp)
set(src ${src} src/TimeStepping/InitTimeStepping.cpp)
set(src ${src} src/TimeStepping/RunTimeStepping.cpp)
set(src ${src} src/TimeStepping/BackwardEuler.cpp)
set(src ${src} src/TimeStepping/CrankNicolson.cpp)
set(src ${src} src/TimeStepping/BackwardEulerAD.cpp)
set(src ${src} src/TimeStepping/CrankNicolsonAD.cpp)


#############################################################
### For math utils                                        ###
#############################################################
set(inc ${inc} include/Utils/MathUtils.h)

### For tensor calculation
set(inc ${inc} include/Utils/RankTwoTensor.h include/Utils/RankFourTensor.h)
set(src ${src} src/Utils/Tensor/RankTwoTensor.cpp)
set(src ${src} src/Utils/Tensor/RankFourTensor.cpp)
#############################################################
### For output system                                     ###
#############################################################
set(inc ${inc} include/OutputSystem/OutputSystem.h)
set(src ${src} src/OutputSystem/OutputSystem.cpp)
set(src ${src} src/OutputSystem/InitOutputStream.cpp)
set(src ${src} src/OutputSystem/WriteResultToVTU.cpp)

#############################################################
### For system level subroutines                          ###
#############################################################
set(inc ${inc} include/FESystem/FESystem.h)
set(src ${src} src/FESystem/InitFESystem.cpp)
set(src ${src} src/FESystem/FESystem.cpp src/FESystem/PrintFESystemInfo.cpp)
set(src ${src} src/FESystem/FormFE.cpp)
set(src ${src} src/FESystem/AssembleLocalToGlobal.cpp)
set(src ${src} src/FESystem/AssembleLocalHistToGlobal.cpp)
set(src ${src} src/FESystem/Projection.cpp)


#############################################################
### For FEProblem(the highest level) subroutines          ###
#############################################################
### Currently, this is the highest level of AsFem
set(inc ${inc} include/FEProblem/FEProblem.h include/FEProblem/FECtrlInfo.h include/FEProblem/FEJobType.h)
set(src ${src} src/FEProblem/RunFEProblem.cpp)
set(src ${src} src/FEProblem/FEProblem.cpp src/FEProblem/ReadInputFile.cpp)
set(src ${src} src/FEProblem/CheckInputFileInfo.cpp)
set(src ${src} src/FEProblem/PreRunFEProblem.cpp)
set(src ${src} src/FEProblem/InitFEProblem.cpp)
set(src ${src} src/FEProblem/StartJob.cpp src/FEProblem/Finalize.cpp)
set(src ${src} src/FEProblem/RunFEAnalysis.cpp)
set(src ${src} src/FEProblem/RunStaticAnalysis.cpp)
set(src ${src} src/FEProblem/RunTransientAnalysis.cpp)
set(src ${src} src/FEProblem/PrintJobInfo.cpp)


##################################################
add_executable(asfem ${inc} ${src})

