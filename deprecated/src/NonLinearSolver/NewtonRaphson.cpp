#include "NonLinearSolver/NonLinearSolver.h"

bool NonLinearSolver::NewtonRaphson(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo){

    _CurrentIters=0;
    bool IsConvergent=false;
    

    
    // feSystem.PresetDirichletBC(mesh,dofHandler,bcSystem,solutionSystem._Unew);
    bcSystem.ApplyPreSetDirichletBC(feCtrlInfo._CurrentT,mesh,dofHandler,solutionSystem._Unew);
    
    if(feCtrlInfo._CurrentStep==0){
        solutionSystem._U=solutionSystem._Unew;
    }
    solutionSystem._V.setZero();
    
    feSystem.ResetMaxAMatrixValue();
    
    while(_CurrentIters<_MaxIters && !IsConvergent){

        if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BACKWARDEULER){
            solutionSystem._Utemp=solutionSystem._Unew;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::CRANKNICOLSON){
            solutionSystem._Utemp=(solutionSystem._U+solutionSystem._Unew)*0.5;
        }
        else{
            solutionSystem._Utemp=solutionSystem._Unew;
        }
        
        
        feSystem.FormFE(feCtrlInfo._ISW,
                        feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,mesh,dofHandler,elmtSystem,mateSystem,
                        solutionSystem._Utemp,solutionSystem._V,solutionSystem._Hist,solutionSystem._HistOld,
                        solutionSystem._ProjValue,
                        equationSystem._AMATRIX,equationSystem._RHS);
        
        // cout<<equationSystem._AMATRIX<<endl;
        bcSystem.SetMaxKMatrixValue(feSystem.GetMaxAMatrixValue()*1.0e5);

        bcSystem.ApplyBC(mesh,dofHandler,feSystem.GetFEPtr(),
                         feCtrlInfo._CurrentT,feCtrlInfo._ctan,
                         equationSystem._AMATRIX,equationSystem._RHS,
                         solutionSystem._Unew);

        
        _Rnorm=equationSystem._RHS.norm()/sqrt(feSystem.GetVolume());
        if(_linearSolver.Solve(equationSystem._AMATRIX,equationSystem._RHS,equationSystem._dU)){
            _Enorm=(equationSystem._RHS*equationSystem._dU).norm()/sqrt(feSystem.GetVolume());


            if(_CurrentIters==0){
                _Rnorm0=_Rnorm;
                _Enorm0=_Enorm;
            }

            _CurrentIters+=1;
            solutionSystem._Unew+=equationSystem._dU;
            feSystem.ResetMaxAMatrixValue();

            bcSystem.ApplyPreSetDirichletBC(feCtrlInfo._CurrentT,mesh,dofHandler,solutionSystem._Unew);

            if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BACKWARDEULER||
               feCtrlInfo._TimeSteppingMethod==TimeSteppingType::CRANKNICOLSON){
                solutionSystem._V=(solutionSystem._Unew-solutionSystem._U)*feCtrlInfo._ctan[1];
            }
            else{
                solutionSystem._V=(solutionSystem._Unew-solutionSystem._U)*feCtrlInfo._ctan[1];
            }


            if(_IsDepDebug){
                PrintDepIterationInfo();
            }

            if(ConvergenceCheck()){
                IsConvergent=true;
                break;
            }
        }
        else{
            cout<<"*** Error: solve Ax=F failed(NewtonRaphson/nonlinearSolver)***"<<endl;
            return false;
        }
    }
    if(_CurrentIters>=_MaxIters && !IsConvergent){
        cout<<"*** Error: newton-raphson solve failed!!!                  ***"<<endl;
        PrintDepIterationInfo();
        return false;
    }
    return IsConvergent;
}