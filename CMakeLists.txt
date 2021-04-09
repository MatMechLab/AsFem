cmake_minimum_required(VERSION 3.8)
project(AsFem)

set(CMAKE_CXX_STANDARD 17)

###############################################
### Set your Eigen/PETSc/MPI path here      ###
### The only things to modify is the        ###
### following two lines(PETSC/MPI_DIR)      ###
###############################################
set(PETSC_DIR "/home/by/Programs/petsc/3.14.3")
set(MPI_DIR "/home/by/Programs/openmpi/4.1.0")

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

###############################################
# For Eigen                                 ###
###############################################
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")

###############################################
### For debug mode                          ###
###############################################
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O2 -fopenmp")

###############################################
### For release mode                        ###
###############################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -O3 -march=native -DNDEBUG")


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
set(src ${src} src/InputSystem/ReadMeshBlock.cpp)
set(src ${src} src/InputSystem/ReadDofsBlock.cpp)
set(src ${src} src/InputSystem/ReadElmtBlock.cpp)
set(src ${src} src/InputSystem/ReadMateBlock.cpp)
set(src ${src} src/InputSystem/ReadBCBlock.cpp)
set(src ${src} src/InputSystem/ReadICBlock.cpp)
set(src ${src} src/InputSystem/ReadQPointBlock.cpp)
set(src ${src} src/InputSystem/ReadOutputBlock.cpp)
set(src ${src} src/InputSystem/ReadPostprocessBlock.cpp)
set(src ${src} src/InputSystem/ReadProjectionBlock.cpp)
set(src ${src} src/InputSystem/ReadNonlinearSolverBlock.cpp)
set(src ${src} src/InputSystem/ReadTimeSteppingBlock.cpp)
set(src ${src} src/InputSystem/ReadFEJobBlock.cpp)

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
### For Mechanics element implementation
set(inc ${inc} include/ElmtSystem/MechanicsElmt.h)
set(src ${src} src/ElmtSystem/MechanicsElmt.cpp)

### for laplace operator
set(src ${src} src/ElmtSystem/LaplaceElmt.cpp)
### for body source term
set(src ${src} src/ElmtSystem/BodySourceElmt.cpp)
### for poisson equation
set(src ${src} src/ElmtSystem/PoissonElmt.cpp)
### for poisson equation
set(src ${src} src/ElmtSystem/TimeDerivElmt.cpp)
### for diffusion equation
set(src ${src} src/ElmtSystem/DiffusionElmt.cpp)
### for diffusion equation
set(src ${src} src/ElmtSystem/CahnHilliardElmt.cpp)
### for Miehe's fracture model
set(src ${src} src/ElmtSystem/MieheFractureElmt.cpp)
### for User1 element
set(src ${src} src/ElmtSystem/User1Elmt.cpp)
### for User2 element
set(src ${src} src/ElmtSystem/User2Elmt.cpp)
### for User3 element
set(src ${src} src/ElmtSystem/User3Elmt.cpp)

#############################################################
### For Materials system in AsFem                         ###
#############################################################
set(inc ${inc} include/MateSystem/MateType.h)
set(inc ${inc} include/MateSystem/MateTypeDefine.h)
set(inc ${inc} include/MateSystem/MateBlock.h)
set(inc ${inc} include/MateSystem/MateSystem.h)
set(src ${src} src/MateSystem/MateSystem.cpp)
### for bulk materials
set(inc ${inc} include/MateSystem/BulkMaterialBase.h)
set(inc ${inc} include/MateSystem/MechanicsMaterialBase.h)
# for elastic material
set(inc ${inc} include/MateSystem/LinearElasticMaterial.h)
### For free energy based materials
set(inc ${inc} include/MateSystem/FreeEnergyMaterialBase.h)
#   for double well free energy materials
set(inc ${inc} include/MateSystem/DoubleWellFreeEnergyMaterial.h)
set(src ${src} src/MateSystem/DoubleWellFreeEnergyMaterial.cpp)
set(inc ${inc} include/MateSystem/BulkMateSystem.h)
set(src ${src} src/MateSystem/BulkMateSystem.cpp)
set(src ${src} src/MateSystem/RunBulkMateLibs.cpp)
set(src ${src} src/MateSystem/ConstPoissonMaterial.cpp)
set(src ${src} src/MateSystem/ConstDiffusionMaterial.cpp)
set(src ${src} src/MateSystem/CahnHilliardMaterial.cpp)
set(src ${src} src/MateSystem/LinearElasticMaterial.cpp)
### for Miehe's phase field fracture materials
set(src ${src} src/MateSystem/MieheFractureMaterial.cpp)
### for UMAT
set(src ${src} src/MateSystem/User1Material.cpp)

#############################################################
### For Boundary condition system in AsFem                ###
#############################################################
set(inc ${inc} include/BCSystem/BCType.h)
set(inc ${inc} include/BCSystem/BCBlock.h)
set(inc ${inc} include/BCSystem/BCSystem.h)
set(src ${src} src/BCSystem/BCSystem.cpp)
set(src ${src} src/BCSystem/CheckAppliedBCNameIsValid.cpp)
set(src ${src} src/BCSystem/PrintBCSystemInfo.cpp)
set(src ${src} src/BCSystem/ApplyBC.cpp)
set(src ${src} src/BCSystem/ApplyDirichletBC.cpp)
set(src ${src} src/BCSystem/ApplyNodalDirichletBC.cpp)
set(src ${src} src/BCSystem/ApplyNeumannBC.cpp)
set(src ${src} src/BCSystem/ApplyNodalNeumannBC.cpp)
set(src ${src} src/BCSystem/RunBCLibs.cpp)
set(src ${src} src/BCSystem/User1BC.cpp)

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
set(src ${src} src/FE/ShapeFuns/LagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFuns/LagrangeShapeFunCalc.cpp)
set(src ${src} src/FE/ShapeFuns/Compute1DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFuns/Compute2DLagrangeShapeFun.cpp)
set(src ${src} src/FE/ShapeFuns/Compute3DLagrangeShapeFun.cpp)
set(inc ${inc} include/FE/ShapeFun.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun.cpp)
#############################################################
### For solution system in AsFem                          ###
#############################################################
set(inc ${inc} include/SolutionSystem/SolutionSystem.h)
set(src ${src} src/SolutionSystem/SolutionSystem.cpp)
set(src ${src} src/SolutionSystem/InitSolution.cpp)

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

#############################################################
### For time stepping system in AsFem                     ###
#############################################################
set(inc ${inc} include/TimeStepping/TimeSteppingType.h)
set(inc ${inc} include/TimeStepping/TimeSteppingBlock.h)
set(inc ${inc} include/TimeStepping/TimeStepping.h)
set(src ${src} src/TimeStepping/TimeStepping.cpp)
set(src ${src} src/TimeStepping/TimeSteppingInit.cpp)
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

