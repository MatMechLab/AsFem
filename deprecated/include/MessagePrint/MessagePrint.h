#ifndef ASFEM_MESSAGEPRINT_H
#define ASFEM_MESSAGEPRINT_H

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

// For program exit
void Msg_AsFem_Exit();

// For input message
void Msg_Input_InvalidArgs();
void Msg_InputFileName_Invalid(string filename);
void Msg_InputFile_NoMeshBlock();
void Msg_InputFile_NoDofsBlock();
// for mesh block
void Msg_Input_LineError(int linenum);
void Msg_Input_BlockBracketNotComplete();
void Msg_Input_Invalid1DMeshType();
void Msg_Input_Invalid2DMeshType();
void Msg_Input_Invalid3DMeshType();

// For dof block
void Msg_Input_DofNameNotFound();
void Msg_Input_DofNameDuplicate();

// For projection block
void Msg_Input_ProjectionNameNotFound();

// For linear solver block
void Msg_Input_LinearSolverBlock_InvalidSolverOption();
void Msg_Input_LinearSolverBlock_MaxIterNotFound();
void Msg_Input_LinearSolverBlock_MaxIterInvalid();
void Msg_Input_LinearSolverBlock_TolNotFound();
void Msg_Input_LinearSolverBlock_TolInvalid();

// For non-linear solver block
void Msg_Input_NonLinearSolverBlock_InvalidSolverOption();
void Msg_Input_NonLinearSolverBlock_MaxIterNotFound();
void Msg_Input_NonLinearSolverBlock_MaxIterInvalid();
void Msg_Input_NonLinearSolverBlock_RrelTolNotFound();
void Msg_Input_NonLinearSolverBlock_RrelTolInvalid();
void Msg_Input_NonLinearSolverBlock_RabsTolNotFound();
void Msg_Input_NonLinearSolverBlock_RabsTolInvalid();
void Msg_Input_NonLinearSolverBlock_ErelTolNotFound();
void Msg_Input_NonLinearSolverBlock_ErelTolInvalid();
void Msg_Input_NonLinearSolverBlock_EabsTolNotFound();
void Msg_Input_NonLinearSolverBlock_EabsTolInvalid();
void Msg_Input_NonLinearSolverBlock_STolNotFound();
void Msg_Input_NonLinearSolverBlock_STolInvalid();

// For qp block
void Msg_Input_QpBlockTypeNotFound();
void Msg_Input_QpBlockOrderNotFound();
void Msg_Input_QpBlockOrderInvalid();

// For BC block
void Msg_Input_BCSubBlockNameNotFound();
void Msg_Input_BCUserNumNotFound();
void Msg_Input_BCUserNumNotValid();
void Msg_Input_BCBoundaryNameNotFound();
void Msg_Input_BCValueNotFound();
void Msg_Input_BCBlockInfoNotComplete();
void Msg_Input_BCBlockNoDof();
void Msg_Input_BCBlockNoElmt();
void Msg_Input_BCBlockNoBoundary();
void Msg_Input_BCBlockUnsupportedType();

void Msg_BCBlockBlockNameInvalid(string blockname);

void Msg_Input_NoBCBlockWarning();

// For IC block
void Msg_Input_ICSubBlockNameNotFound();
void Msg_Input_ICUserNumNotFound();
void Msg_Input_ICUserNumNotValid();
void Msg_Input_ICBlockNameNotFound();
void Msg_Input_ICValueNotFound();
void Msg_Input_ICBlockInfoNotComplete();
void Msg_Input_ICBlockNoDof();
void Msg_Input_ICBlockNoElmt();
void Msg_Input_ICBlockNoBlock();
void Msg_Input_ICBlockNoValueWarning(string blockname);
void Msg_Input_ICBlockUnsupportedType();
void Msg_ICBlockBlockNameInvalid(string blockname);


// For Elmt block
void Msg_Input_ElmtUserNumNotFound();
void Msg_Input_ElmtUserNumNotValid();
void Msg_Input_ElmtSubBlockNameNotFound();
void Msg_Input_ElmtBlockUnsupportedType();
void Msg_Input_ElmtBlockMateNameNotFound();
void Msg_Input_ElmtBlockBlockNameNotFound();
void Msg_Input_ElmtBlockDuplicateDofsName();



void Msg_Input_ElmtBlockInfoNotComplete();
void Msg_Input_ElmtBlockNoType();
void Msg_Input_ElmtBlockNoDofs();

void Msg_ElmtBlockBlockNameInvalid(string blockname);

void Msg_InputFile_NoElmtBlockFound();
void Msg_InputFile_NoICBlockFoundWarning();

// For material block
void Msg_Input_MateSubBlockNameNotFound();
void Msg_Input_MateUserNumNotFound();
void Msg_Input_MateUserNumNotValid();
void Msg_Input_MateParamsNotFound();
void Msg_Input_MateBlockInfoNotComplete();
void Msg_Input_MateBlockNoElmt();
void Msg_Input_MateBlockNoParams();
void Msg_Input_MateBlockNoParamsWarning(string blockname);
void Msg_Input_MateBlockUnsupportType();

void Msg_InputFile_NoMateBlockFoundWarning();





void Msg_Input_NoDofNameFound();

// For Job block
void Msg_Input_JobBlock_InvalidJobType();
void Msg_Input_JobBlock_InvalidDebugOption();
void Msg_Input_JobBlock_InvalidAdaptiveOption();
void Msg_Input_JobBlock_InvalidProjectionOption();
void Msg_Input_JobBlock_InvalidSolverOption();
void Msg_Input_JobBlock_InvalidNonLinearSolverOption();
void Msg_Input_JobBlock_InvalidTimeSteppingOption();
void Msg_Input_JobBlock_NoTimeFound();
void Msg_Input_JobBlock_NoDtFound();
void Msg_Input_JobBlock_InvalidTimeValue();
void Msg_Input_JobBlock_InvalidDtValue();
void Msg_Input_JobBblock_TransientInfoNotComplete();
void Msg_Input_JobBlock_NoIntervalFound();
void Msg_Input_JobBlock_InvalidIntervalValue();
void Msg_InputFile_NoJobBlockFound();
//************************************************
//*** Message for input file check
void Msg_InputCheck_DofInvalidInElmtBlock(int i,string blockname);
void Msg_InputCheck_DofInvalidInBCBlock(string blockname);
void Msg_InputCheck_DofInvalidInICBlock(string blockname);
void Msg_InputCheck_MateNotFoundForElmtBlock(string matename,string blockname);
void Msg_InputCheck_ElmtBlockDuplicateDomain(string str);
void Msg_InputCheck_DofIsNotAssignedToElmt(string dofname);
//*****************************************************
//*** Message for gauss point generation
//*****************************************************
void Msg_GaussPoint_Invalid1DLintNum(int n);
void Msg_GaussPoint_Invalid1DLobbotLintNum(int n);

void Msg_GaussPoint_Invalid2DGaussOrder(int n);
void Msg_GaussPoint_Invalid2DMeshType();
void Msg_GaussPoint_Invalid2DLobattoMeshType();

void Msg_GaussPoint_Invalid3DGaussOrder(int n);
void Msg_GaussPoint_Invalid3DMeshType();
void Msg_GaussPoint_Invalid3DLobattoMeshType();

//***********************************************
//**** Message for shape function print       ***
//***********************************************
void Msg_ShapeFun_InformationNotComlete();
void Msg_ShapeFun_UnsupportedMeshType();

void Msg_ShapeFun_Lagrange1DElmtSingular();
void Msg_ShapeFun_Lagrange2DElmtSingular();
void Msg_ShapeFun_Lagrange3DElmtSingular();


//***********************************************
//**** Message for linear solver print        ***
//***********************************************
void Msg_LinearSolver_AnalyzeFailed();
void Msg_LinearSolver_FactorizeFailed();
void Msg_LinearSolver_SolveFailed();

//***********************************************
//**** Message for equation system print      ***
//***********************************************
void Msg_EquationSystem_NotInitialized();

//***********************************************
//**** Message for nonlinear solver print     ***
//***********************************************
void Msg_NonLinearSolver_UnsupportSolver();


//***********************************************
//**** Message for material system print      ***
//***********************************************
void Msg_MaterialSystem_UnsupportUmat();

#endif