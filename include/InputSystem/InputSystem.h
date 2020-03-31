//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_INPUTSYSTEM_H
#define ASFEM_INPUTSYSTEM_H


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>


#include "petsc.h"

//**************************************
//*** For AsFem's own header file
//**************************************
#include "Utils/StringUtils.h"
#include "MessagePrinter/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"
#include "Solution/Solution.h"
#include "FE/FE.h"

#include "FEProblem/JobBlock.h"

#include "NonlinearSolver/NonlinearSolverBlock.h"
#include "TimeStepping/TimeSteppingBlock.h"

#include "OutputSystem/OutputBlock.h"

class InputSystem{
public:
    InputSystem(int args,char *argv[]);
    InputSystem();
    void InitInputSystem(int args,char *argv[]);
    bool ReadInputFile(Mesh &mesh,
                       DofHandler &dofHandler,
                       ElmtSystem &elmtSystem,
                       MateSystem &mateSystem,
                       BCSystem &bcSystem,
                       ICSystem &icSystem,
                       Solution &solution,
                       FE &fe,
                       NonlinearSolverBlock &nonlinearSolverBlock,
                       TimeSteppingBlock &timesteppingblock,
                       JobBlock &jobBlock,
                       OutputBlock &outputblock);

private:
    //**************************************
    //*** functions for reading each block
    //**************************************
    bool ReadMeshBlock(ifstream &in,string str,int &linenum,Mesh &mesh);
    bool ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler);
    bool ReadElmtBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler);
    bool ReadMateBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MateSystem &mateSystem);

    bool ReadQPointBlock(ifstream &in,string str,int &linenum,FE &fe);
    //*** for boundary blocks
    bool ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem,DofHandler &dofHandler);

    //*** for initial condition block
    bool ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem,DofHandler &dofHandler);

    //*** for job block
    bool ReadJobBlock(ifstream &in,string str,int &linenum,JobBlock &jobBlock);

    //*** for nonlinear solver block
    bool ReadNonlinearSolverBlock(ifstream &in,string str,int &linenum,NonlinearSolverBlock &nonlinearSolverBlock);

    //*** for time stepping block
    bool ReadTimeSteppingBlock(ifstream &in,string str,int &linenum,TimeSteppingBlock &timesteppingblock);

    //*** for projection block
    bool ReadProjectionBlock(ifstream &in,string str,int &linenum,Solution &solution);

    //*** for output block
    bool ReadOutputBlock(ifstream &in,string str,int &linenum,OutputBlock &outputblock);

public:
    //**************************************
    //*** Basic getting functioins
    //**************************************
    inline string GetInputFileName() const {return _InputFileName;}
    inline string GetMeshFileName() const {return _MeshFileName;}

private:
    string _InputFileName,_MeshFileName;
    bool _IsBuiltInMesh=true;
    bool _HasInputFileName;
};


#endif //ASFEM_INPUTSYSTEM_H