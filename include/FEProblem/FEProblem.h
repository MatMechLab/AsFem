//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_FEPROBLEM_H
#define ASFEM_FEPROBLEM_H


#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <chrono>

#include "petsc.h"

//***********************************
//*** For AsFem's own header file
//***********************************
#include "MessagePrinter/MessagePrinter.h"

#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"

#include "Solution/Solution.h"
#include "EquationSystem/EquationSystem.h"

#include "FE/FE.h"
#include "FESystem/FESystem.h"

#include "JobBlock.h"

#include "NonlinearSolver/NonlinearSolver.h"

#include "TimeStepping/TimeStepping.h"

#include "OutputSystem/OutputSystem.h"
#include "FeCtrlInfo.h"

class FEProblem{
public:
    FEProblem(int args,char *argv[]);

    void RunFEProblem();

private:
    bool RunInputSystem();
    void InitFEProblem();
    void RunAnalysis();

private:
    //********************************************
    //*** For static and transient analysis
    //********************************************
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

    Solution _solution;
    EquationSystem _equationSystem;

    FE _fe;
    FESystem _feSystem;

    JobBlock _jobBlock;

    NonlinearSolver _nonlinearsolver;

    TimeStepping _timestepping;

    OutputSystem _outputSystem;

    TimeSteppingBlock _timesteppingblock;
    NonlinearSolverBlock _nonlinearsolverblock;

    FeCtrlInfo _feCtrlInfo;

private:
    //****************************************
    //*** for time duration of each module
    //****************************************
    PetscMPIInt _rank;
    PetscReal Duration(chrono::high_resolution_clock::time_point &p1,chrono::high_resolution_clock::time_point &p2){
        return chrono::duration_cast<std::chrono::microseconds>(p2-p1).count()/1.0e6;
    }

    // timer count for each module
    chrono::high_resolution_clock::time_point _TimerStartOfInput,_TimerEndOfInput;
    PetscReal _DurationOfInput;
    chrono::high_resolution_clock::time_point _TimerStartOfInputCheck,_TimerEndOfInputChekc;
    PetscReal _DurationOfInputCheck;

    chrono::high_resolution_clock::time_point _TimerStartOfElmtInit,_TimerEndOfElmtInit;
    PetscReal _DurationOfElmtInit;

    chrono::high_resolution_clock::time_point _TimerStartOfDofInit,_TimerEndOfDofInit;
    PetscReal _DurationOfDofInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFormFE,_TimerEndOfFormFE;
    PetscReal _DurationOfFormFE;
    chrono::high_resolution_clock::time_point _TimerStartOfOutput,_TimerEndOfOutput;
    PetscReal _DurationOfOutput;
    chrono::high_resolution_clock::time_point _TimerStartOfSparseInit,_TimerEndOfSparseInit;
    PetscReal _DurationOfSparseInit;

    chrono::high_resolution_clock::time_point _TimerStartOfSolInit,_TimerEndOfSolInit;
    PetscReal _DurationOfSolInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFESysInit,_TimerEndOfFESysInit;
    PetscReal _DurationOfFESysInit;

    chrono::high_resolution_clock::time_point _TimerStartOfSolverInit,_TimerEndOfSolverInit;
    PetscReal _DurationOfSolverInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFEProInit,_TimerEndOfFEProInit;
    PetscReal _DurationOfFEProInit;

    

    chrono::system_clock::time_point _JobStartTime,_JobEndTime;
    string _JobStartTimeContext,_JobEndTimeContext;

    chrono::high_resolution_clock::time_point _TimerStartOfJob,_TimerEndOfJob;
    double _DurationOfJob;

};


#endif // ASFEM_FEPROBLEM_H