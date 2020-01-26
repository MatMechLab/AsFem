#include "NonLinearSolver/NonLinearSolver.h"

NonLinearSolver::NonLinearSolver(NonLinearSolverType method){
    if(method==NonLinearSolverType::LINEAR){
        _SolverType=NonLinearSolverType::LINEAR;
    }
    else if(method==NonLinearSolverType::NEWTONRAPHSON){
        _SolverType=NonLinearSolverType::NEWTONRAPHSON;
    }
    else if(method==NonLinearSolverType::MODIFIEDNEWTON){
        _SolverType=NonLinearSolverType::MODIFIEDNEWTON;
    }
    else if(method==NonLinearSolverType::NEWTONWITHLINESEARCH){
        _SolverType=NonLinearSolverType::NEWTONWITHLINESEARCH;
    }
    else if(method==NonLinearSolverType::BFGS){
        _SolverType=NonLinearSolverType::BFGS;
    }
    else{
        Msg_NonLinearSolver_UnsupportSolver();
        _SolverType=NonLinearSolverType::NEWTONRAPHSON;
    }

    _IsDebug=false;
    _ConvergenceCriterion=0;//0-->use the residual or energy
                            //1-->only use the residual
                            //2-->only use the energy
                            //3-->use both the residual and energy

    _MaxIters=100;_CurrentIters=0;
    _RAbsTol=1.0e-8;_RRelTol=1.0e-9;
    _EAbsTol=1.0e-19;_ERelTol=1.0e-21;
    _Rnorm0=1.0;_Rnorm=1.0;
    _Enorm0=1.0;_Enorm=1.0;
    _Stol=1.0e-4;
    _MaxLineSearch=500;
}