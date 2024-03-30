#############################################################
#############################################################
### For beginners, please don't edit the following line!  ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
### Do not edit the following lines !!!                   ###
#############################################################
#############################################################
# For main.cpp
set(inc "")
set(src src/main.cpp)

#############################################################
### For message printer utils                             ###
#############################################################
set(inc ${inc} include/Utils/MessagePrinter.h include/Utils/MessageColor.h)
set(src ${src} src/Utils/MessagePrinter.cpp)

#############################################################
### For string utils                                      ###
#############################################################
set(inc ${inc} include/Utils/StringUtils.h)
set(src ${src} src/Utils/StringUtils.cpp)

#############################################################
### For json utils                                        ###
#############################################################
set(inc ${inc} include/Utils/JsonUtils.h)
set(src ${src} src/Utils/JsonUtils.cpp)

#############################################################
### For timer                                             ###
#############################################################
set(inc ${inc} include/Utils/Timer.h)
set(src ${src} src/Utils/Timer.cpp)

#############################################################
### For MPITool                                           ###
#############################################################
set(inc ${inc} include/MPIUtils/MPIDataBus.h)
set(src ${src} src/MPIUtils/MPIDataBus.cpp)

#############################################################
### For mathematic utils                                  ###
#############################################################
set(inc ${inc} include/MathUtils/Vector3d.h)
set(src ${src} src/MathUtils/Vector3d.cpp)
set(inc ${inc} include/MathUtils/VectorXd.h)
set(src ${src} src/MathUtils/VectorXd.cpp)
set(inc ${inc} include/MathUtils/MatrixXd.h)
set(src ${src} src/MathUtils/MatrixXd.cpp)
### for PETSc vector
set(inc ${inc} include/MathUtils/Vector.h)
set(src ${src} src/MathUtils/Vector.cpp)
### for PETSc's sparse matrix
set(inc ${inc} include/MathUtils/SparseMatrix.h)
set(src ${src} src/MathUtils/SparseMatrix.cpp)
### for rank-2 tensor
set(inc ${inc} include/MathUtils/Rank2Tensor.h)
set(src ${src} src/MathUtils/Rank2Tensor.cpp)
### for rank-4 tensor
set(inc ${inc} include/MathUtils/Rank4Tensor.h)
set(src ${src} src/MathUtils/Rank4Tensor.cpp)
### for general math funs
set(inc ${inc} include/MathUtils/MathFuns.h)
set(src ${src} src/MathUtils/MathFuns.cpp)


#############################################################
### For mesh class                                        ###
#############################################################
### for nodes
set(inc ${inc} include/Mesh/Nodes.h)
### for mesh type
set(inc ${inc} include/Mesh/MeshType.h)
### for mesh structure data
set(inc ${inc} include/Mesh/MeshData.h)
### for mesh generator class
set(inc ${inc} include/Mesh/MeshGeneratorBase.h)
### for 1d mesh generator class
set(inc ${inc} include/Mesh/Lagrange1DMeshGenerator.h)
set(src ${src} src/Mesh/Lagrange1DMeshGenerator.cpp)
### for 2d mesh generator class
set(inc ${inc} include/Mesh/Lagrange2DMeshGenerator.h)
set(src ${src} src/Mesh/Lagrange2DMeshGenerator.cpp)
### for 3d mesh generator class
set(inc ${inc} include/Mesh/Lagrange3DMeshGenerator.h)
set(src ${src} src/Mesh/Lagrange3DMeshGenerator.cpp)
### for general mesh generator class
set(inc ${inc} include/Mesh/MeshGenerator.h)
set(src ${src} src/Mesh/MeshGenerator.cpp)
### for bulk mesh
set(inc ${inc} include/Mesh/BulkMesh.h)
set(src ${src} src/Mesh/BulkMesh.cpp)
set(src ${src} src/Mesh/SaveMesh.cpp)
set(src ${src} src/Mesh/PrintMesh.cpp)
### for the final mesh
set(inc ${inc} include/Mesh/Mesh.h)
### utils for msh file
set(inc ${inic} include/Mesh/MshFileUtils.h)
set(src ${src} src/Mesh/MshFileUtils.cpp)
### for msh2 file
set(inc ${inc} include/Mesh/Msh2FileImporter.h)
set(src ${src} src/Mesh/Msh2FileImporter.cpp)
### for msh2 file
set(inc ${inc} include/Mesh/Msh4FileImporter.h)
set(src ${src} src/Mesh/Msh4FileImporter.cpp)
### for gmsh2 file
set(inc ${inc} include/Mesh/Gmsh2FileImporter.h)
set(src ${src} src/Mesh/Gmsh2FileImporter.cpp)
### for MeshIO
set(inc ${inc} include/Mesh/MeshFileImporterBase.h)
set(inc ${inc} include/Mesh/MeshFileImporter.h)
set(src ${src} src/Mesh/MeshFileImporter.cpp)

#############################################################
### For FE cell class                                     ###
#############################################################
### for single FE mesh cell
set(inc ${inc} include/FECell/SingleMeshCell.h)
set(inc ${inc} include/FECell/FECellData.h)
set(inc ${inc} include/FECell/FECell.h)
set(src ${src} src/FECell/FECell.cpp)
set(inc ${inc} include/FECell/FECellGeneratorBase.h)
### for 1d lagrange mesh cell
### for edge2
set(inc ${inc} include/FECell/Lagrange1DEdge2MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange1DEdge2MeshCellGenerator.cpp)
### for edge3
set(inc ${inc} include/FECell/Lagrange1DEdge3MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange1DEdge3MeshCellGenerator.cpp)
### for edge4
set(inc ${inc} include/FECell/Lagrange1DEdge4MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange1DEdge4MeshCellGenerator.cpp)
### for 2d lagrange mesh cell
### for quad4
set(inc ${inc} include/FECell/Lagrange2DQuad4MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange2DQuad4MeshCellGenerator.cpp)
### for quad8
set(inc ${inc} include/FECell/Lagrange2DQuad8MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange2DQuad8MeshCellGenerator.cpp)
### for quad9
set(inc ${inc} include/FECell/Lagrange2DQuad9MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange2DQuad9MeshCellGenerator.cpp)
### for 3d lagrange mesh cell
### for 3d-hex8 mesh cell
set(inc ${inc} include/FECell/Lagrange3DHex8MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange3DHex8MeshCellGenerator.cpp)
### for 3d-hex20 mesh cell
set(inc ${inc} include/FECell/Lagrange3DHex20MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange3DHex20MeshCellGenerator.cpp)
### for 3d-hex27 mesh cell
set(inc ${inc} include/FECell/Lagrange3DHex27MeshCellGenerator.h)
set(src ${src} src/FECell/Lagrange3DHex27MeshCellGenerator.cpp)
##
set(inc ${inc} include/FECell/FECellGenerator.h)
set(src ${src} src/FECell/FECellGenerator.cpp)

#############################################################
### For FE cell importer class                            ###
#############################################################
set(inc ${inc} include/FECell/MeshFile2FECellImporterBase.h)
### for msh2 file importer
set(inc ${inc} include/FECell/Msh2File2FECellImporter.h)
set(src ${src} src/FECell/Msh2File2FECellImporter.cpp)


#############################################################
### For inputystem                                        ###
#############################################################
set(inc ${inc} include/InputSystem/InputSystem.h)
set(src ${src} src/InputSystem/InputSystem.cpp)
set(src ${src} src/InputSystem/ReadInputFile.cpp)
### for mesh block
set(src ${src} src/InputSystem/ReadMeshBlock.cpp)
### for dofs block
set(src ${src} src/InputSystem/ReadDofsBlock.cpp)
### for dofs block
set(src ${src} src/InputSystem/ReadElmtsBlock.cpp)
### for qpoint block
set(src ${src} src/InputSystem/ReadQPointBlock.cpp)
### for shapefunction block
set(src ${src} src/InputSystem/ReadShapeFunBlock.cpp)
### for boundary condition block
set(src ${src} src/InputSystem/ReadBCsBlock.cpp)
### for initial condition block
set(src ${src} src/InputSystem/ReadICsBlock.cpp)
### for projection condition block
set(src ${src} src/InputSystem/ReadProjectionBlock.cpp)
### for nonlinear solver block
set(src ${src} src/InputSystem/ReadNLSolverBlock.cpp)
### for nonlinear solver block
set(src ${src} src/InputSystem/ReadTimeSteppingBlock.cpp)
### for output block
set(src ${src} src/InputSystem/ReadOutputBlock.cpp)
### for postprocess block
set(src ${src} src/InputSystem/ReadPostprocessBlock.cpp)
### for job block
set(src ${src} src/InputSystem/ReadJobBlock.cpp)


#############################################################
### For dof handlers                                      ###
#############################################################
set(inc ${inc} include/DofHandler/DofHandler.h)
### for bulk dofHandler
set(inc ${inc} include/DofHandler/BulkDofHandler.h)
set(src ${src} src/DofHandler/BulkDofHandler.cpp)
set(src ${src} src/DofHandler/BulkDofHandlerSettings.cpp)
set(src ${src} src/DofHandler/CreateBulkDofsMap.cpp)

#############################################################
### For boundary conditions                               ###
#############################################################
set(inc ${inc} include/BCSystem/BCType.h)
set(inc ${inc} include/BCSystem/BCBlock.h)
set(inc ${inc} include/BCSystem/BCSystem.h)
set(src ${src} src/BCSystem/BCSystem.cpp)
set(src ${src} src/BCSystem/ApplyBoundaryConditions.cpp)
set(src ${src} src/BCSystem/RunBCLibs.cpp)
set(src ${src} src/BCSystem/BCSystemAssemble.cpp)
### for dirichlet bcs
set(inc ${inc} include/BCSystem/DirichletBCBase.h)
###
set(inc ${inc} include/BCSystem/DirichletBC.h)
set(src ${src} src/BCSystem/DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/RotatedDirichletBC.h)
set(src ${src} src/BCSystem/RotatedDirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/User1DirichletBC.h)
set(src ${src} src/BCSystem/User1DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/User2DirichletBC.h)
set(src ${src} src/BCSystem/User2DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/User3DirichletBC.h)
set(src ${src} src/BCSystem/User3DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/User4DirichletBC.h)
set(src ${src} src/BCSystem/User4DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/User5DirichletBC.h)
set(src ${src} src/BCSystem/User5DirichletBC.cpp)
###
set(inc ${inc} include/BCSystem/Poisson2DBenchmarkBC.h)
set(src ${src} src/BCSystem/Poisson2DBenchmarkBC.cpp)
###
set(src ${src} src/BCSystem/ApplyDirichletBC.cpp)
### for integrated bcs
set(src ${src} src/BCSystem/ApplyIntegratedBC.cpp)
set(inc ${inc} include/BCSystem/IntegrateBCBase.h)
###
set(inc ${inc} include/BCSystem/NeumannBC.h)
set(src ${src} src/BCSystem/NeumannBC.cpp)
### for pressure bc
set(inc ${inc} include/BCSystem/PressureBC.h)
set(src ${src} src/BCSystem/PressureBC.cpp)
### for traction bc
set(inc ${inc} include/BCSystem/TractionBC.h)
set(src ${src} src/BCSystem/TractionBC.cpp)


#############################################################
### For initial conditions                                ###
#############################################################
set(inc ${inc} include/ICSystem/ICType.h)
set(inc ${inc} include/ICSystem/ICBlock.h)
### ic base class
set(inc ${inc} include/ICSystem/InitialConditionBase.h)
### for constant ic
set(inc ${inc} include/ICSystem/ConstantIC.h)
set(src ${src} src/ICSystem/ConstantIC.cpp)
###
set(inc ${inc} include/ICSystem/RandomIC.h)
set(src ${src} src/ICSystem/RandomIC.cpp)
###
set(inc ${inc} include/ICSystem/CircleIC.h)
set(src ${src} src/ICSystem/CircleIC.cpp)
###
set(inc ${inc} include/ICSystem/ICSystem.h)
set(src ${src} src/ICSystem/ICSystem.cpp)
set(src ${src} src/ICSystem/ApplyInitialConditions.cpp)
set(src ${src} src/ICSystem/RunICLibs.cpp)

#############################################################
### For element system                                    ###
#############################################################
set(inc ${inc} include/ElmtSystem/ElmtType.h)
set(inc ${inc} include/ElmtSystem/ElmtBlock.h)
set(inc ${inc} include/ElmtSystem/LocalElmtData.h)
### for bulk element system
set(inc ${inc} include/ElmtSystem/BulkElmtSystem.h)
set(src ${src} src/ElmtSystem/BulkElmtSystem.cpp)
set(src ${src} src/ElmtSystem/BulkElmtSystemInit.cpp)
set(src ${src} src/ElmtSystem/RunBulkElmtLibs.cpp)
set(inc ${inc} include/ElmtSystem/BulkElmtBase.h)
### for poisson element
set(inc ${inc} include/ElmtSystem/PoissonElement.h)
set(src ${src} src/ElmtSystem/PoissonElement.cpp)
### for diffusion element
set(inc ${inc} include/ElmtSystem/DiffusionElement.h)
set(src ${src} src/ElmtSystem/DiffusionElement.cpp)
###
set(inc ${inc} include/ElmtSystem/AllenCahnElement.h)
set(src ${src} src/ElmtSystem/AllenCahnElement.cpp)
###
set(inc ${inc} include/ElmtSystem/MechanicsElement.h)
set(src ${src} src/ElmtSystem/MechanicsElement.cpp)
###
set(inc ${inc} include/ElmtSystem/CahnHilliardElement.h)
set(src ${src} src/ElmtSystem/CahnHilliardElement.cpp)
###
set(inc ${inc} include/ElmtSystem/KobayashiElement.h)
set(src ${src} src/ElmtSystem/KobayashiElement.cpp)
###
set(inc ${inc} include/ElmtSystem/StressDiffusionElement.h)
set(src ${src} src/ElmtSystem/StressDiffusionElement.cpp)
###
set(inc ${inc} include/ElmtSystem/AllenCahnFractureElement.h)
set(src ${src} src/ElmtSystem/AllenCahnFractureElement.cpp)
###
set(inc ${inc} include/ElmtSystem/MieheFractureElement.h)
set(src ${src} src/ElmtSystem/MieheFractureElement.cpp)
###
set(inc ${inc} include/ElmtSystem/StressCahnHilliardElement.h)
set(src ${src} src/ElmtSystem/StressCahnHilliardElement.cpp)
### for scalar body source element
set(inc ${inc} include/ElmtSystem/ScalarBodySourceElement.h)
set(src ${src} src/ElmtSystem/ScalarBodySourceElement.cpp)
### for laplacian element
set(inc ${inc} include/ElmtSystem/LaplaceElement.h)
set(src ${src} src/ElmtSystem/LaplaceElement.cpp)
### for diffusion allencahn fracture
set(inc ${inc} include/ElmtSystem/DiffusionACFractureElement.h)
set(src ${src} src/ElmtSystem/DiffusionACFractureElement.cpp)
### for element system
set(inc ${inc} include/ElmtSystem/ElmtSystem.h)
set(src ${src} src/ElmtSystem/ElmtSystem.cpp)

#############################################################
### For material system                                   ###
#############################################################
set(inc ${inc} include/MateSystem/MateType.h)
set(inc ${inc} include/MateSystem/MaterialsName.h)
### for material container
set(inc ${inc} include/MateSystem/MaterialsContainer.h)
set(src ${src} src/MateSystem/MaterialsContainer.cpp)
### for bulk material system
set(inc ${inc} include/MateSystem/BulkMateSystem.h)
set(inc ${inc} include/MateSystem/BulkMaterialBase.h)
set(src ${src} src/MateSystem/BulkMateSystem.cpp)
set(src ${src} src/MateSystem/InitBulkMateLibs.cpp)
set(src ${src} src/MateSystem/RunBulkMateLibs.cpp)
### for poisson material
set(inc ${inc} include/MateSystem/ConstPoissonMaterial.h)
set(src ${src} src/MateSystem/ConstPoissonMaterial.cpp)
### for poisson 1d+2d benchmark material
set(inc ${inc} include/MateSystem/Poisson1DBenchmarkMaterial.h)
set(src ${src} src/MateSystem/Poisson1DBenchmarkMaterial.cpp)
set(inc ${inc} include/MateSystem/Poisson2DBenchmarkMaterial.h)
set(src ${src} src/MateSystem/Poisson2DBenchmarkMaterial.cpp)
### for nonlinear poisson2d material
set(inc ${inc} include/MateSystem/NonlinearPoisson2DMaterial.h)
set(src ${src} src/MateSystem/NonlinearPoisson2DMaterial.cpp)
### for nonlinear poisson3d material
set(inc ${inc} include/MateSystem/NonlinearPoisson3DMaterial.h)
set(src ${src} src/MateSystem/NonlinearPoisson3DMaterial.cpp)
### for const diffusion material
set(inc ${inc} include/MateSystem/ConstDiffusionMaterial.h)
set(src ${src} src/MateSystem/ConstDiffusionMaterial.cpp)
### for const diffusion material
set(inc ${inc} include/MateSystem/NonlinearDiffusion2DMaterial.h)
set(src ${src} src/MateSystem/NonlinearDiffusion2DMaterial.cpp)
###
set(inc ${inc} include/MateSystem/DoubleWellPotentialMaterial.h)
set(src ${src} src/MateSystem/DoubleWellPotentialMaterial.cpp)
### for elastic/plastic material
set(inc ${inc} include/MateSystem/ElasticMaterialBase.h)
set(inc ${inc} include/MateSystem/PlasticMaterialBase.h)
### for linear elastic material
set(inc ${inc} include/MateSystem/LinearElasticMaterial.h)
set(src ${src} src/MateSystem/LinearElasticMaterial.cpp)
### for saint venant hyperelastic material
set(inc ${inc} include/MateSystem/SaintVenantMaterial.h)
set(src ${src} src/MateSystem/SaintVenantMaterial.cpp)
### for saint venant hyperelastic material
set(inc ${inc} include/MateSystem/NeoHookeanMaterial.h)
set(src ${src} src/MateSystem/NeoHookeanMaterial.cpp)
### for elastoplastic material
set(inc ${inc} include/MateSystem/SmallStrainJ2PlasticityMaterial.h)
set(src ${src} src/MateSystem/SmallStrainJ2PlasticityMaterial.cpp)
###
set(inc ${inc} include/MateSystem/SmallStrainExpLawJ2PlasticityMaterial.h)
set(src ${src} src/MateSystem/SmallStrainExpLawJ2PlasticityMaterial.cpp)
### for free energy material
set(inc ${inc} include/MateSystem/FreeEnergyMaterialBase.h)
###
set(inc ${inc} include/MateSystem/BinaryMixtureMaterial.h)
set(src ${src} src/MateSystem/BinaryMixtureMaterial.cpp)
### kobayasi dendrite material
set(inc ${inc} include/MateSystem/KobayashiDendriteMaterial.h)
set(src ${src} src/MateSystem/KobayashiDendriteMaterial.cpp)
### for coupled materials
set(inc ${inc} include/MateSystem/SmallStrainDiffusionMaterial.h)
set(src ${src} src/MateSystem/SmallStrainDiffusionMaterial.cpp)
###
set(inc ${inc} include/MateSystem/LinearElasticFractureMaterial.h)
set(src ${src} src/MateSystem/LinearElasticFractureMaterial.cpp)
###
set(inc ${inc} include/MateSystem/MieheFractureMaterial.h)
set(src ${src} src/MateSystem/MieheFractureMaterial.cpp)
###
set(inc ${inc} include/MateSystem/NeoHookeanPFFractureMaterial.h)
set(src ${src} src/MateSystem/NeoHookeanPFFractureMaterial.cpp)
###
set(inc ${inc} include/MateSystem/SmallStrainCahnHilliardMaterial.h)
set(src ${src} src/MateSystem/SmallStrainCahnHilliardMaterial.cpp)
###
set(inc ${inc} include/MateSystem/SmallStrainDiffusionJ2Material.h)
set(src ${src} src/MateSystem/SmallStrainDiffusionJ2Material.cpp)
###
set(inc ${inc} include/MateSystem/DiffusionACFractureMaterial.h)
set(src ${src} src/MateSystem/DiffusionACFractureMaterial.cpp)
### User-1 material
set(inc ${inc} include/MateSystem/User1Material.h)
set(src ${src} src/MateSystem/User1Material.cpp)
### for material system
set(inc ${inc} include/MateSystem/MateSystem.h)
set(src ${src} src/MateSystem/MateSystem.cpp)

#############################################################
### For FE space class                                    ###
#############################################################
set(inc ${inc} include/FE/FE.h)
set(src ${src} src/FE/FE.cpp)
####################################
# for shape function data
####################################
set(inc ${inc} include/FE/ShapeFunType.h)
### for 1d shape function base 
set(inc ${inc} include/FE/ShapeFun1DBase.h)
### for 1d shape function edge2 
set(inc ${inc} include/FE/ShapeFun1DEdge2.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun1DEdge2.cpp)
### for 1d shape function edge3 
set(inc ${inc} include/FE/ShapeFun1DEdge3.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun1DEdge3.cpp)
### for 1d shape function edge4 
set(inc ${inc} include/FE/ShapeFun1DEdge4.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun1DEdge4.cpp)
### for 1d shape functions
set(inc ${inc} include/FE/ShapeFun1D.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun1D.cpp)
#####################################
### for 2d shape function  
#####################################
set(inc ${inc} include/FE/ShapeFun2DBase.h)
### for 2d tri3 shape function
set(inc ${inc} include/FE/ShapeFun2DTri3.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2DTri3.cpp)
### for 2d tri6 shape function
set(inc ${inc} include/FE/ShapeFun2DTri6.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2DTri6.cpp)
### for 2d quad4 shape function
set(inc ${inc} include/FE/ShapeFun2DQuad4.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2DQuad4.cpp)
### for 2d quad8 shape function
set(inc ${inc} include/FE/ShapeFun2DQuad8.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2DQuad8.cpp)
### for 2d quad9 shape function
set(inc ${inc} include/FE/ShapeFun2DQuad9.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2DQuad9.cpp)
### for 2d shape function
set(inc ${inc} include/FE/ShapeFun2D.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun2D.cpp)
#####################################
### for 3d shape function 
#####################################
set(inc ${inc} include/FE/ShapeFun3DBase.h)
### for 3d tet4 shape function
set(inc ${inc} include/FE/ShapeFun3DTet4.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3DTet4.cpp)
### for 3d tet10 shape function
set(inc ${inc} include/FE/ShapeFun3DTet10.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3DTet10.cpp)
### for 3d hex8 shape function
set(inc ${inc} include/FE/ShapeFun3DHex8.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3DHex8.cpp)
### for 3d hex20 shape function
set(inc ${inc} include/FE/ShapeFun3DHex20.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3DHex20.cpp)
### for 3d hex27 shape function
set(inc ${inc} include/FE/ShapeFun3DHex27.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3DHex27.cpp)
### for 3d shape function
set(inc ${inc} include/FE/ShapeFun3D.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun3D.cpp)
#####################################
### for user-defined shape functions 
#####################################
set(inc ${inc} include/FE/ShapeFunUserBase.h)
### user1 shape function
set(inc ${inc} include/FE/ShapeFunUser1.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser1.cpp)
### user2 shape function
set(inc ${inc} include/FE/ShapeFunUser2.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser2.cpp)
### user3 shape function
set(inc ${inc} include/FE/ShapeFunUser3.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser3.cpp)
### user4 shape function
set(inc ${inc} include/FE/ShapeFunUser4.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser4.cpp)
### user5 shape function
set(inc ${inc} include/FE/ShapeFunUser5.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser5.cpp)
### user shape function
set(inc ${inc} include/FE/ShapeFunUser.h)
set(src ${src} src/FE/ShapeFuns/ShapeFunUser.cpp)
### for general shape functions
set(inc ${inc} include/FE/ShapeFun.h)
set(src ${src} src/FE/ShapeFuns/ShapeFun.cpp)
set(src ${src} src/FE/ShapeFuns/ShapeFunCalc.cpp)
####################################
# for QPoint classes
####################################
set(inc ${inc} include/FE/QPointType.h)
set(inc ${inc} include/FE/QPointGeneratorBase.h)
### for 1d qpoint generator
set(inc ${inc} include/FE/QPoint1DGenerator.h)
set(src ${src} src/FE/QPoints/QPoint1DGenerator.cpp)
### for 2d qpoint generator
set(inc ${inc} include/FE/QPoint2DGenerator.h)
set(src ${src} src/FE/QPoints/QPoint2DGenerator.cpp)
### for 3d qpoint generator
set(inc ${inc} include/FE/QPoint3DGenerator.h)
set(src ${src} src/FE/QPoints/QPoint3DGenerator.cpp)
### for 1d lobatto generator
set(inc ${inc} include/FE/QPoint1DLobattoGenerator.h)
set(src ${src} src/FE/QPoints/QPoint1DLobattoGenerator.cpp)
### for 2d lobatto generator
set(inc ${inc} include/FE/QPoint2DLobattoGenerator.h)
set(src ${src} src/FE/QPoints/QPoint2DLobattoGenerator.cpp)
### for 3d lobatto generator
set(inc ${inc} include/FE/QPoint3DLobattoGenerator.h)
set(src ${src} src/FE/QPoints/QPoint3DLobattoGenerator.cpp)
### for user-defined qpoint generation
set(inc ${inc} include/FE/QPointUserGeneratorBase.h)
### for user-1 qpoints generator
set(inc ${inc} include/FE/QPointUser1Generator.h)
set(src ${src} src/FE/QPoints/QPointUser1Generator.cpp)
### for user-2 qpoints generator
set(inc ${inc} include/FE/QPointUser2Generator.h)
set(src ${src} src/FE/QPoints/QPointUser2Generator.cpp)
### for user-3 qpoints generator
set(inc ${inc} include/FE/QPointUser3Generator.h)
set(src ${src} src/FE/QPoints/QPointUser3Generator.cpp)
### for user-4 qpoints generator
set(inc ${inc} include/FE/QPointUser4Generator.h)
set(src ${src} src/FE/QPoints/QPointUser4Generator.cpp)
### for user-5 qpoints generator
set(inc ${inc} include/FE/QPointUser5Generator.h)
set(src ${src} src/FE/QPoints/QPointUser5Generator.cpp)
### for final gauss point class
set(inc ${inc} include/FE/QPoint.h)
set(src ${src} src/FE/QPoints/QPoint.cpp)

#############################################################
### For FESystem class                                    ###
#############################################################
set(inc ${inc} include/FESystem/FECalcType.h)
set(inc ${inc} include/FESystem/FESystem.h)
### for bulk FE system
set(inc ${inc} include/FESystem/BulkFESystem.h)
set(src ${src} src/FESystem/BulkFESystem.cpp)
set(src ${src} src/FESystem/BulkFESystemInit.cpp)
set(src ${src} src/FESystem/FormBulkFE.cpp)
set(src ${src} src/FESystem/BulkFESystemAssemble.cpp)
### for FE system
set(src ${src} src/FESystem/FESystem.cpp)

#############################################################
### For Projection system                                 ###
#############################################################
set(inc ${inc} include/ProjectionSystem/ProjectionType.h)
set(inc ${inc} include/ProjectionSystem/ProjectionData.h)
set(inc ${inc} include/ProjectionSystem/ProjectionBase.h)
set(inc ${inc} include/ProjectionSystem/ProjectionSystem.h)
### for leastsquare projection
set(inc ${inc} include/ProjectionSystem/LeastSquareProjection.h)
set(src ${src} src/ProjectionSystem/LeastSquareProjection.cpp)
###
set(src ${src} src/ProjectionSystem/ProjectionSystem.cpp)
set(src ${src} src/ProjectionSystem/Projection.cpp)
set(src ${src} src/ProjectionSystem/ProjectionLibs.cpp)
set(src ${src} src/ProjectionSystem/ProjectionGettings.cpp)

#############################################################
### For Equation class                                    ###
#############################################################
set(inc ${inc} include/EquationSystem/EquationSystem.h)
set(src ${src} src/EquationSystem/EquationSystem.cpp)

#############################################################
### For Solution class                                    ###
#############################################################
set(inc ${inc} include/SolutionSystem/SolutionSystem.h)
set(src ${src} src/SolutionSystem/SolutionSystem.cpp)
set(src ${src} src/SolutionSystem/SolutionSystemInit.cpp)
set(src ${src} src/SolutionSystem/SolutionUpdate.cpp)

#############################################################
### For Nonlinear solver class                            ###
#############################################################
set(inc ${inc} include/NonlinearSolver/NonlinearSolverBase.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolverType.h)
set(inc ${inc} include/NonlinearSolver/NonlinearSolverBlock.h)
### for SNES solver
set(inc ${inc} include/NonlinearSolver/SNESSolver.h)
set(src ${src} src/NonlinearSolver/SNESSolver.cpp)
set(src ${src} src/NonlinearSolver/SNESSolve.cpp)
###
set(inc ${inc} include/NonlinearSolver/NonlinearSolver.h)
set(src ${src} src/NonlinearSolver/NonlinearSolver.cpp)

#############################################################
### For Timesteppinig class                               ###
#############################################################
set(inc ${inc} include/TimeStepping/TimeSteppingType.h)
set(inc ${inc} include/TimeStepping/TimeSteppingData.h)
set(inc ${inc} include/TimeStepping/TimeSteppingTool.h)
###
set(inc ${inc} include/TimeStepping/TimeStepping.h)
set(src ${src} src/TimeStepping/TimeStepping.cpp)
set(src ${src} src/TimeStepping/TimeSteppingSolve.cpp)

#############################################################
### For Result output class                               ###
#############################################################
set(inc ${inc} include/OutputSystem/ResultWriterBase.h)
set(inc ${inc} include/OutputSystem/ResultFileFormat.h)
set(inc ${inc} include/OutputSystem/OutputSystem.h)
set(src ${src} src/OutputSystem/OutputSystem.cpp)
set(src ${src} src/OutputSystem/SaveResults.cpp)
### for vtu format
set(inc ${inc} include/OutputSystem/VTUWriter.h)
set(src ${src} src/OutputSystem/VTUWriter.cpp)


#############################################################
### For Postprocessor class                               ###
#############################################################
set(inc ${inc} include/Postprocess/Postprocessor.h)
set(inc ${inc} include/Postprocess/PostprocessorBlock.h)
set(inc ${inc} include/Postprocess/PostprocessorType.h)
set(src ${src} src/Postprocess/Postprocessor.cpp)
set(src ${src} src/Postprocess/SavePostprocessResults.cpp)
set(src ${src} src/Postprocess/ExecutePostprocess.cpp)
###
set(inc ${inc} include/Postprocess/NodalPostprocessorBase.h)
###
set(inc ${inc} include/Postprocess/NodalValuePostprocessor.h)
set(src ${src} src/Postprocess/NodalValuePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/NodalScalarMatePostprocessor.h)
set(src ${src} src/Postprocess/NodalScalarMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/NodalVectorMatePostprocessor.h)
set(src ${src} src/Postprocess/NodalVectorMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/NodalRank2MatePostprocessor.h)
set(src ${src} src/Postprocess/NodalRank2MatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/NodalRank4MatePostprocessor.h)
set(src ${src} src/Postprocess/NodalRank4MatePostprocessor.cpp)
###
set(src ${src} src/Postprocess/ExecuteNodalPostprocess.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralPostprocessorBase.h)
set(src ${src} src/Postprocess/ExecuteSideIntegralPostprocess.cpp)
###
set(src ${src} src/Postprocess/RunSideIntegralPostprocessLibs.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralValuePostprocessor.h)
set(src ${src} src/Postprocess/SideIntegralValuePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/AreaPostprocessor.h)
set(src ${src} src/Postprocess/AreaPostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralScalarMatePostprocessor.h)
set(src ${src} src/Postprocess/SideIntegralScalarMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralVectorMatePostprocessor.h)
set(src ${src} src/Postprocess/SideIntegralVectorMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralRank2MatePostprocessor.h)
set(src ${src} src/Postprocess/SideIntegralRank2MatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/SideIntegralRank4MatePostprocessor.h)
set(src ${src} src/Postprocess/SideIntegralRank4MatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/User1SideIntegralPostprocessor.h)
set(src ${src} src/Postprocess/User1SideIntegralPostprocessor.cpp)
### for volume integral pps
set(inc ${inc} include/Postprocess/VolumeIntegralPostprocessorBase.h)
###
set(inc ${inc} include/Postprocess/VolumePostprocessor.h)
set(src ${src} src/Postprocess/VolumePostprocessor.cpp)
###
set(src ${src} src/Postprocess/RunVolumeIntegralPostprocessLibs.cpp)
set(src ${src} src/Postprocess/ExecuteVolumeIntegralPostprocess.cpp)
###
set(inc ${inc} include/Postprocess/VolumeIntegralValuePostprocessor.h)
set(src ${src} src/Postprocess/VolumeIntegralValuePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/VolumeIntegralScalarMatePostprocessor.h)
set(src ${src} src/Postprocess/VolumeIntegralScalarMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/VolumeIntegralVectorMatePostprocessor.h)
set(src ${src} src/Postprocess/VolumeIntegralVectorMatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/VolumeIntegralRank2MatePostprocessor.h)
set(src ${src} src/Postprocess/VolumeIntegralRank2MatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/VolumeIntegralRank4MatePostprocessor.h)
set(src ${src} src/Postprocess/VolumeIntegralRank4MatePostprocessor.cpp)
###
set(inc ${inc} include/Postprocess/User1VolumeIntegralPostprocessor.h)
set(src ${src} src/Postprocess/User1VolumeIntegralPostprocessor.cpp)

#############################################################
### For FEProblem class                                   ###
#############################################################
set(inc ${inc} include/FEProblem/FEProblem.h)
set(inc ${inc} include/FEProblem/FEControlInfo.h)
set(inc ${inc} include/FEProblem/FEJobType.h)
set(inc ${inc} include/FEProblem/FEJobBlock.h)
set(src ${src} src/FEProblem/FEProblem.cpp)
set(src ${src} src/FEProblem/RunFEProblem.cpp)
set(src ${src} src/FEProblem/RunStaticAnalysis.cpp)
set(src ${src} src/FEProblem/RunTransientAnalysis.cpp) 


#############################################################
### For Application class                                 ###
#############################################################
set(inc ${inc} include/Application.h)
set(src ${src} src/Application.cpp)