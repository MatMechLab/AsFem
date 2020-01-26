#include "TimeStepping/TimeStepping.h"

void TimeStepping::RunTimeStepping(Mesh &mesh,DofHandler &dofHandler,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            FESystem &feSystem,NonLinearSolver &nonLinearSolver,
            FECtrlInfo &feCtrlInfo,
            OutputSystem &outputSystem){
    _TimerStart=chrono::high_resolution_clock::now();
    if(!_IsAdaptive){
        switch(_Method){
        case TimeSteppingType::BACKWARDEULER:
            BackwardEuler(mesh,dofHandler,solutionSystem,equationSystem,
                          elmtSystem,mateSystem,bcSystem,icSystem,
                          feSystem,nonLinearSolver,feCtrlInfo,outputSystem);
            break;
        case TimeSteppingType::CRANKNICOLSON:
            CrankNicolson(mesh,dofHandler,solutionSystem,equationSystem,
                          elmtSystem,mateSystem,bcSystem,icSystem,
                          feSystem,nonLinearSolver,feCtrlInfo,outputSystem);
            break;
        default:
            break;
        }
    }
    else{
        switch(_Method){
        case TimeSteppingType::BACKWARDEULER:
            BackwardEulerAD(mesh,dofHandler,solutionSystem,equationSystem,
                          elmtSystem,mateSystem,bcSystem,icSystem,
                          feSystem,nonLinearSolver,feCtrlInfo,outputSystem);
            break;
        case TimeSteppingType::CRANKNICOLSON:
            CrankNicolsonAD(mesh,dofHandler,solutionSystem,equationSystem,
                          elmtSystem,mateSystem,bcSystem,icSystem,
                          feSystem,nonLinearSolver,feCtrlInfo,outputSystem);
            break;
        default:
            break;
        }
    }
}