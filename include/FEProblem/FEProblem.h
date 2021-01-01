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
//+++ Date   : 2020.06.29
//+++ Purpose: Define the FE problem analysis class in AsFem,
//+++          It is designed as the top level of the whole AsFem
//+++          framework
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <ctime>
#include <chrono>

#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "FE/FE.h"
#include "FESystem/FESystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"
#include "NonlinearSolver/NonlinearSolver.h"
#include "TimeStepping/TimeStepping.h"
#include "OutputSystem/OutputSystem.h"

#include "FEProblem/FEJobType.h"
#include "FEProblem/FEControlInfo.h"
#include "FEProblem/FEJobBlock.h"

using namespace std;

class FEProblem{
public:
    FEProblem();

    void InitFEProblem(int args,char *argv[]);

    void Run();

    void Finalize();

private:
    void ReadInputFile();
    void InitAllComponents();

    void RunStaticAnalysis();
    void RunTransientAnalysis();

private:
    InputSystem _inputSystem;
    Mesh _mesh;
    DofHandler _dofHandler;
    ElmtSystem _elmtSystem;
    MateSystem _mateSystem;
    BCSystem _bcSystem;
    ICSystem _icSystem;
    FE _fe;
    FESystem _feSystem;
    SolutionSystem _solutionSystem;
    EquationSystem _equationSystem;
    NonlinearSolver _nonlinearSolver;
    TimeStepping _timestepping;
    OutputSystem _outputSystem;

    FEJobType _feJobType;

    FEControlInfo _feCtrlInfo;

    FEJobBlock _feJobBlock;

private:
    //****************************************************************
    //*** for profiling
    //****************************************************************
    PetscMPIInt _rank;
    double Duration(chrono::high_resolution_clock::time_point &p1,chrono::high_resolution_clock::time_point &p2){
        return chrono::duration_cast<std::chrono::microseconds>(p2-p1).count()/1.0e6;
    }

    chrono::high_resolution_clock::time_point _TimerStart,_TimerEnd;
    double _Duration;


};