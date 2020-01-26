#include "NonLinearSolver/NonLinearSolver.h"

bool NonLinearSolver::NewtonRaphsonWithLineSearch(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo){

    _CurrentIters=0;
    bool IsConvergent=false;
    bool IsLineSearchFound=false;
    

    
    // feSystem.PresetDirichletBC(mesh,dofHandler,bcSystem,solutionSystem._Unew);
    bcSystem.ApplyPreSetDirichletBC(feCtrlInfo._CurrentT,mesh,dofHandler,solutionSystem._Unew);
    if(feCtrlInfo._CurrentStep==0){
        solutionSystem._U=solutionSystem._Unew;
    }
    solutionSystem._V.setZero();
    
    _Rnorm=0.0;
    
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

        bcSystem.SetMaxKMatrixValue(feSystem.GetMaxAMatrixValue()*1.0e3);

        bcSystem.ApplyBC(mesh,dofHandler,feSystem.GetFEPtr(),
                         feCtrlInfo._CurrentT,feCtrlInfo._ctan,
                         equationSystem._AMATRIX,equationSystem._RHS,
                         solutionSystem._Unew);
        // cout<<equationSystem._AMATRIX<<endl;
        _Rnorm=equationSystem._RHS.norm()/sqrt(feSystem.GetVolume());
        if(_linearSolver.Solve(equationSystem._AMATRIX,equationSystem._RHS,equationSystem._dU)){
            _Enorm=(equationSystem._RHS*equationSystem._dU).norm()/sqrt(feSystem.GetVolume());
            _s0=AtDotB(equationSystem._RHS,equationSystem._dU)/sqrt(feSystem.GetVolume());
            _eta0=1.0;
            _s=_s0;

            if(_CurrentIters==0){
                _Rnorm0=_Rnorm;
                _Enorm0=_Enorm;
            }
            // start to do line search to find the suitable update
            IsLineSearchFound=false;_nLineSearch=0;
            _eta=_eta0;
            while(_s/_s0>_Stol && _nLineSearch<_MaxLineSearch){
                solutionSystem._Utemp=solutionSystem._Unew+equationSystem._dU*_eta0;
                // only forms the residual
                feSystem.FormFE(3,
                        feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,mesh,dofHandler,elmtSystem,mateSystem,
                        solutionSystem._Utemp,solutionSystem._V,solutionSystem._Hist,solutionSystem._HistOld,
                        solutionSystem._ProjValue,
                        equationSystem._AMATRIX,equationSystem._RHS);
                bcSystem.ApplyBC(mesh,dofHandler,feSystem.GetFEPtr(),
                         feCtrlInfo._CurrentT,feCtrlInfo._ctan,
                         equationSystem._AMATRIX,equationSystem._RHS,
                         solutionSystem._Utemp);
                // calculate new s
                _s=AtDotB(equationSystem._RHS,equationSystem._dU)/sqrt(feSystem.GetVolume());
                // then we can calculate the new eta
                _eta=_eta0*_s0/(_s0-_s);
                _eta0=_eta;
                _nLineSearch+=1;

                if(_IsDepDebug){
                    printf("***    Iters=%3d, %2d search: s=%11.3e,eta=%11.3e ***\n",_CurrentIters+1,_nLineSearch,_s,_eta);
                }

                if(_s/_s0<=_Stol){
                    IsLineSearchFound=true;
                    break;
                }
            }
            if(IsLineSearchFound){
                solutionSystem._Unew+=equationSystem._dU*_eta;
                bcSystem.ApplyPreSetDirichletBC(feCtrlInfo._CurrentT,mesh,dofHandler,solutionSystem._Unew);
                if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BACKWARDEULER||
                   feCtrlInfo._TimeSteppingMethod==TimeSteppingType::CRANKNICOLSON){
                    solutionSystem._V=(solutionSystem._Unew-solutionSystem._U)*feCtrlInfo._ctan[1];
                }
                else{
                    solutionSystem._V=(solutionSystem._Unew-solutionSystem._U)*feCtrlInfo._ctan[1];
                }

                feSystem.ResetMaxAMatrixValue();
                _CurrentIters+=1;
            }
            else{
                cout<<"*** Error: line search failed in NR iteration !!!          ***"<<endl;
                PrintDepIterationInfo();
                return false;
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