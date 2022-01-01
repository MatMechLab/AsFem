//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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
#include "InputSystem/DofsBlockReader.h"
#include "InputSystem/ElmtsBlockReader.h"
#include "InputSystem/MatesBlockReader.h"
#include "InputSystem/BCsBlockReader.h"
#include "InputSystem/ICsBlockReader.h"
#include "InputSystem/QPointBlockReader.h"
#include "InputSystem/OutputBlockReader.h"
#include "InputSystem/PostprocessBlockReader.h"
#include "InputSystem/ProjectionBlockReader.h"
#include "InputSystem/NonlinearSolverBlockReader.h"
#include "InputSystem/TimeSteppingBlockReader.h"
#include "InputSystem/FEJobBlockReader.h"


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


class InputSystem:public MeshBlockReader,
                  public DofsBlockReader,
                  public ElmtsBlockReader,
                  public MatesBlockReader,
                  public BCsBlockReader,
                  public ICsBlockReader,
                  public QPointBlockReader,
                  public OutputBlockReader,
                  public PostprocessBlockReader,
                  public ProjectionBlockReader,
                  public NonlinearSolverBlockReader,
                  public TimeSteppingBlockReader,
                  public FEJobBlockReader
{
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
    //*** private variables
    //******************************************************
    MeshIO _meshio;
    NonlinearSolverBlock _nonlinearSolverBlock;
    string _InputFileName,_MeshFileName;
    bool _HasInputFileName=false;
    bool _IsBuiltInMesh=true;
    bool _IsReadOnly=false;

};
