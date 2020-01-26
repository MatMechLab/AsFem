#ifndef ASFEM_TIMESTEPPING_H
#define ASFEM_TIMESTEPPING_H

#include <iostream>
#include <iomanip>
#include <chrono>

#include "TimeSteppingType.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCs/BCSystem.h"
#include "ICs/ICSystem.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MaterialSystem/MaterialSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "FESystem/FESystem.h"
#include "FEProblem/FECtrlInfo.h"
#include "NonLinearSolver/NonLinearSolver.h"
#include "OutputSystem/OutputSystem.h"

using namespace std;

class TimeStepping
{
public:
    TimeStepping();

    void InitTimeStepping();

    void SetAdaptiveFlag(bool flag){_IsAdaptive=flag;}
    void SetTimeSteppingMethod(TimeSteppingType method){_Method=method;}
    void SetGrowthFactor(const double &factor) {_GrowthFactor=factor;}
    void SetCutbackFactor(const double &factor) {_CutbackFactor=factor;}
    void SetMaxDt(const double &dt) {_MaxDt=dt;}
    void SetMinDt(const double &dt) {_MinDt=dt;}

    inline double GetGrowthFactor() const {return _GrowthFactor;}
    inline double GetCutbackFactor() const {return _CutbackFactor;}
    void SetOptimIters(const int iters){_OptIters=iters;}
    inline int GetOptimIters() const{return _OptIters;}
    inline double GetMaxDt() const {return _MaxDt;}

    double ctan[2];

    void RunTimeStepping(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem);
private:
    void AdaptiveStepping();
    void ConstStepping();

// time integration algorithm
private:
    void BackwardEuler(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem);
    // adaptive time stepping
    void BackwardEulerAD(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem);
    void CrankNicolson(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem);
    void CrankNicolsonAD(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem);

private:
    bool _IsAdaptive=false;
    TimeSteppingType _Method;
    double _TotalTime,_dt0,_dt;

    // for time adaptive
    double _GrowthFactor=1.2,_CutbackFactor=0.5;
    int _OptIters=3;
    double _MaxDt,_MinDt;

    long int _TotalStep,_CurrentStep;
    double _CurrentTime=0.0;

    chrono::high_resolution_clock::time_point _TimerStart,_TimerEnd;
    double _Duration;

};

#endif // ASFEM_TIMESTEPPING_H