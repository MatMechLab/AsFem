#ifndef ASFEM_INPUTSYSTEM_H
#define ASFEM_INPUTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>


// AsFem's own header file
#include "MessagePrint/MessagePrint.h"
#include "Utils/StringUtils.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MaterialSystem/MaterialSystem.h"
#include "FE/QPBlockInfo.h"
#include "FEProblem/FECtrlInfo.h"

#include "LinearSolver/LinearSolverBlockInfo.h"
#include "NonLinearSolver/NonLinearSolverBlockInfo.h"

#include "SolutionSystem/SolutionSystem.h"

using namespace std;


class InputSystem
{
public:
    InputSystem(int args,char *argv[]);
    InputSystem();
    inline string GetInputFileName() const {return _InputFileName;}
    void InitInputSystem(int args,char *argv[]);
    bool ReadInputFile(Mesh &mesh,
                       DofHandler &dofHandler,
                       BCSystem &bcSystem,
                       ICSystem &icSystem,
                       ElmtSystem &elmtSystem,
                       MaterialSystem &mateSystem,
                       QPBlockInfo &qpBlockInfo,
                       LinearSolverBlockInfo &linearSolverBlockInfo,
                       NonLinearSolverBlockInfo &nonlinearSolverBlockInfo,
                       SolutionSystem &solutionSystem,
                       FECtrlInfo &feCtrlInfo);

private:
    bool ReadMeshBlock(ifstream &in,string str,int &linenum,Mesh &mesh);
    bool ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler);

    bool ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem);
    
    bool ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem);

    bool ReadElmtBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler);
    bool ReadMateBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MaterialSystem &mateSystem);

    bool ReadQpBlock(ifstream &in,string str,int &linenum,QPBlockInfo &qpBlockInfo);
    
    // For linear solver block
    bool ReadLinearSolverBlock(ifstream &in,string str,int &linenum,LinearSolverBlockInfo &linearSolverBlockInfo);
    // For nonlinear solver block
    bool ReadNonLinearSolverBlock(ifstream &in,string str,int &linenum,NonLinearSolverBlockInfo &nonlinearSolverBlockInfo);


    // For projection block
    bool ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem);

    bool ReadJobBlock(ifstream &in,string str,int &linenum,FECtrlInfo &feCtrlInfo);
    

private:
    string _InputFileName;
    bool _HasInputFileName=false;

    // information for mesh block
    string _MeshFileName="";
    bool _UseBultiInMesh=true;

};

#endif 