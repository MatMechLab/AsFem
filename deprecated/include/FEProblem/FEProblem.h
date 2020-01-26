#ifndef ASFEM_FEPROBLEM_H
#define ASFEM_FEPROBLEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <chrono>

// For AsFem's built-in classes
#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MaterialSystem/MaterialSystem.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "LinearSolver/LinearSolverBlockInfo.h"
#include "LinearSolver/LinearSolver.h"

#include "NonLinearSolver/NonLinearSolverBlockInfo.h"
#include "NonLinearSolver/NonLinearSolver.h"

#include "FESystem/FESystem.h"
#include "TimeStepping/TimeStepping.h"
#include "OutputSystem/OutputSystem.h"
#include "FEJobType.h"
#include "FECtrlInfo.h"
#include "FE/QPBlockInfo.h"


using namespace std;


class FEProblem
{
public:
    FEProblem(int args,char *argv[]);

    void RunFEProblem();

private:
    void StartJob();
    void PreRunFEProblem();
    void ReadInputFile();
    void CheckInputFileInfo();
    //******************************
    //*** create dof maps and generate sparsity patterns
    void InitFEProblem();
    //******************************
    void RunFEAnalysis();
    void RunStaticAnalysis();
    void RunTransientAnalysis();
    void Finalize();

    void PrintJobInfo();


private:
    chrono::system_clock::time_point _JobStartTime,_JobEndTime;
    string _JobStartTimeContext,_JobEndTimeContext;
    // timer count for each module
    chrono::high_resolution_clock::time_point _TimerStartOfInput,_TimerEndOfInput;
    double _DurationOfInput;
    chrono::high_resolution_clock::time_point _TimerStartOfInputCheck,_TimerEndOfInputChekc;
    double _DurationOfInputCheck;

    chrono::high_resolution_clock::time_point _TimerStartOfElmtInit,_TimerEndOfElmtInit;
    double _DurationOfElmtInit;

    chrono::high_resolution_clock::time_point _TimerStartOfDofInit,_TimerEndOfDofInit;
    double _DurationOfDofInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFormFE,_TimerEndOfFormFE;
    double _DurationOfFormFE;
    chrono::high_resolution_clock::time_point _TimerStartOfOutput,_TimerEndOfOutput;
    double _DurationOfOutput;
    chrono::high_resolution_clock::time_point _TimerStartOfSparseInit,_TimerEndOfSparseInit;
    double _DurationOfSparseInit;

    chrono::high_resolution_clock::time_point _TimerStartOfSolInit,_TimerEndOfSolInit;
    double _DurationOfSolInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFESysInit,_TimerEndOfFESysInit;
    double _DurationOfFESysInit;

    chrono::high_resolution_clock::time_point _TimerStartOfSolverInit,_TimerEndOfSolverInit;
    double _DurationOfSolverInit;

    chrono::high_resolution_clock::time_point _TimerStartOfFEProInit,_TimerEndOfFEProInit;
    double _DurationOfFEProInit;

    chrono::high_resolution_clock::time_point _TimerStartOfStaticJob,_TimerEndOfStaticJob;
    double _DurationOfStaticJob;

private:
    InputSystem inputSystem;
    Mesh mesh;
    BCSystem bcSystem;
    ICSystem icSystem;
    DofHandler dofHandler;
    ElmtSystem elmtSystem;
    MaterialSystem materialSystem;
    QPBlockInfo qpBlockInfo;

    EquationSystem equationSystem;
    SolutionSystem solutionSystem;

    FESystem feSystem;

    TimeStepping timeStepping;

    LinearSolver linearSolver;
    NonLinearSolver nonlinearSolver;

    LinearSolverBlockInfo linearSolverBlock;
    NonLinearSolverBlockInfo nonlinearSolverBlock;

    FEJobType JobType=FEJobType::STATIC;

    FECtrlInfo feCtrlInfo;
    OutputSystem outputSystem;

};

#endif // ASFEM_FEPROBLEM_H