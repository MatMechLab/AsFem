//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Implement the input file reader for AsFem
//+++          So, in this class, the input file for AsFem is only
//+++          accessible via this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>


#include "petsc.h"

/**
 * For the utils of message print and string tools
 */
#include "Utils/StringUtils.h"
#include "Utils/MessagePrinter.h"

/**
 * For different block readers related to each block in the input file
 */
#include "InputSystem/MeshBlockReader.h"

#include "Mesh/Mesh.h"
#include "Mesh/MeshIO.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"
#include "FE/FE.h"
#include "SolutionSystem/SolutionSystem.h"
#include "NonlinearSolver/NonlinearSolver.h"
#include "TimeStepping/TimeStepping.h"
#include "OutputSystem/OutputSystem.h"
#include "Postprocess/Postprocess.h"
#include "FEProblem/FEJobBlock.h"


class InputSystem:public MeshBlockReader{
public:
    InputSystem(int args,char *argv[]);
    InputSystem();
    void InitInputSystem(int args,char *argv[]);

    bool ReadInputFile(Mesh &mesh,DofHandler &dofHandler,ElmtSystem &elmtSystem,MateSystem &mateSystem,
                       BCSystem &bcSystem,ICSystem &icSystem,
                       FE &fe,
                       SolutionSystem &solutionSystem,
                       OutputSystem &outputSystem,
                       Postprocess &postProcessSystem,
                       NonlinearSolver &nonlinearSolver,
                       TimeStepping &timestepping,
                       FEJobBlock &feJobBlock);

    bool IsReadOnlyMode()const{return _IsReadOnly;}

private:
    //******************************************************
    //*** functions for reading each block
    //******************************************************
    bool ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler);
    bool ReadElmtBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler);
    bool ReadMateBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MateSystem &mateSystem);

    //******************************************************
    //*** functions for reading bcs and ics
    //******************************************************
    bool ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem,DofHandler &dofHandler);
    bool ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem,DofHandler &dofHandler);

    //******************************************************
    //*** functions for reading [qpoint]
    //******************************************************
    bool ReadQPointBlock(ifstream &in,string str,int &linenum,FE &fe);

    //******************************************************
    //*** functions for reading [output]
    //******************************************************
    bool ReadOutputBlock(ifstream &in,string str,int &linenum,OutputSystem &outputSystem);

    //******************************************************
    //*** functions for reading [postprocess]
    //******************************************************
    bool ReadPostprocessBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,Postprocess &postprocess,DofHandler &dofHandler);

    //******************************************************
    //*** functions for reading [projection]
    //******************************************************
    bool ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem);

    //******************************************************
    //*** functions for reading [nonlinearsolver]
    //******************************************************
    bool ReadNonlinearSolverBlock(ifstream &in,string str,int &linenum,NonlinearSolver &nonlinearSolver);


    //******************************************************
    //*** functions for reading [timestepping]
    //******************************************************
    bool ReadTimeSteppingBlock(ifstream &in,string str,int &linenum,TimeStepping &timestepping);

    //******************************************************
    //*** functions for reading [job]
    //******************************************************
    bool ReadFEJobBlock(ifstream &in,string str,int &linenum,FEJobBlock &feJobBlock);

    
    //******************************************************
    //*** private variables
    //******************************************************
    MeshIO _meshio;
    NonlinearSolverBlock _nonlinearSolverBlock;
    string _InputFileName,_MeshFileName;
    bool _HasInputFileName=false;
    bool _IsBuiltInMesh=true;
    bool _IsReadOnly=false;

};
