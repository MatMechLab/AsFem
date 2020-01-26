#ifndef ASFEM_FECTRLINFO_H
#define ASFEM_FECTRLINFO_H

#include <iomanip>
#include "FEJobType.h"
#include "LinearSolver/LinearSolverType.h"
#include "NonLinearSolver/NonLinearSolverType.h"
#include "TimeStepping/TimeSteppingType.h"

class FECtrlInfo
{
public:
    bool _IsDebugOn=false;
    bool _IsDebugDep=false;
    double _CurrentT=0.0;
    double _CurrentDt=0.0;
    long int _CurrentStep=0;
    bool _IsProjOn=false;
    int _OutputInterval=1;
    FEJobType _JobType=FEJobType::STATIC;

    bool _IsPrintMesh=false;

    TimeSteppingType _TimeSteppingMethod=TimeSteppingType::BACKWARDEULER;
    int _ISW=6;
    int _Interval=1;
    bool _IsAdaptive=false;
    double _dt0=1.0e-6;
    double _dt=1.0e-6,dtold=1.0e-6;
    double _EndTime=1.0e10;
    double _DtMax=0.1;
    double _DtMin=1.0e-10;
    int _OptIters=3;
    double _ctan[2]={1.0,1.0};
    
};

#endif // ASFEM_FECTRLINFO_H