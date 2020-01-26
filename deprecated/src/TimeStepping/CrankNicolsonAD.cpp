#include "TimeStepping/TimeStepping.h"

void TimeStepping::CrankNicolsonAD(Mesh &mesh,DofHandler &dofHandler,
        SolutionSystem &solutionSystem,EquationSystem &equationSystem,
        ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
        BCSystem &bcSystem,ICSystem &icSystem,
        FESystem &feSystem,NonLinearSolver &nonLinearSolver,
        FECtrlInfo &feCtrlInfo,
        OutputSystem &outputSystem){
    // init the basic settings
    bool IsConvergent=false;
    _dt0=feCtrlInfo._dt0;
    _dt=feCtrlInfo._dt0;
    _CurrentTime=0.0;
    _TotalTime=feCtrlInfo._EndTime;

    feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BACKWARDEULER;
    feCtrlInfo._ISW=6;

    feCtrlInfo._ctan[0]=0.5;
    feCtrlInfo._ctan[1]=1.0/feCtrlInfo._dt0;

    icSystem.ApplyIC(mesh,dofHandler,solutionSystem._U);
    solutionSystem._Unew=solutionSystem._U;

    // initial history variable
    feSystem.FormFE(4,0.0,0.0,feCtrlInfo._ctan,
                    mesh,dofHandler,
                    elmtSystem,mateSystem,
                    solutionSystem._Unew,solutionSystem._V,
                    solutionSystem._Hist,solutionSystem._HistOld,
                    solutionSystem._ProjValue,
                    equationSystem._AMATRIX,equationSystem._RHS);
    solutionSystem._HistOld=solutionSystem._Hist;


    _TotalStep=(int)(feCtrlInfo._EndTime/feCtrlInfo._dt0);
    feCtrlInfo._dt=feCtrlInfo._dt0;
    _CurrentStep=0;
    if(feCtrlInfo._IsProjOn){
        outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U,
                                      solutionSystem.GetProjNumPerNode(),
                                      solutionSystem.GetProjNameList(),
                                      solutionSystem._ProjValue);
    }
    else{
        outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U);
    }
    printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());

    for(_CurrentTime=0.0;_CurrentTime<=_TotalTime;){
        feCtrlInfo._ctan[1]=1.0/feCtrlInfo._dt;
        _CurrentStep+=1;
        _CurrentTime+=feCtrlInfo._dt;
        feCtrlInfo._CurrentStep=_CurrentStep;
        feCtrlInfo._CurrentT=_CurrentTime;
        if(feCtrlInfo._IsDebugOn || feCtrlInfo._IsDebugDep){
                printf("*** Step=%7d, t=%12.6e, dt=%12.6e          ***\n",(int)_CurrentStep,_CurrentTime,feCtrlInfo._dt);
        }
        if(nonLinearSolver.Solve(mesh,dofHandler,solutionSystem,equationSystem,
                                 elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo)){
            solutionSystem._U=solutionSystem._Unew;
            solutionSystem.UpdateSolution(_CurrentStep);
            // update history variables
            feSystem.FormFE(8,feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,
                            mesh,dofHandler,elmtSystem,mateSystem,
                            solutionSystem._U,solutionSystem._V,
                            solutionSystem._Hist,solutionSystem._HistOld,
                            solutionSystem._ProjValue,
                            equationSystem._AMATRIX,equationSystem._RHS);
            solutionSystem._HistOld=solutionSystem._Hist;

            if(_CurrentStep%feCtrlInfo._Interval==0){
                if(feCtrlInfo._IsDebugOn && !feCtrlInfo._IsDebugDep){
                    nonLinearSolver.PrintIterationInfo();
                }
                if(feCtrlInfo._IsProjOn){
                    feSystem.FormFE(9,feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,
                            mesh,dofHandler,elmtSystem,mateSystem,
                            solutionSystem._U,solutionSystem._V,
                            solutionSystem._Hist,solutionSystem._HistOld,
                            solutionSystem._ProjValue,
                            equationSystem._AMATRIX,equationSystem._RHS);
                    outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U,
                                          solutionSystem.GetProjNumPerNode(),
                                          solutionSystem.GetProjNameList(),
                                          solutionSystem._ProjValue);
                }
                else{
                    outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U);
                }
                printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());
            }
        }
        else{
            // if current solve failed
            _CurrentStep-=1;
            _CurrentTime-=feCtrlInfo._dt;
            feCtrlInfo._dt*=_CutbackFactor;
            feCtrlInfo._CurrentStep=_CurrentStep;
            feCtrlInfo._CurrentT=_CurrentTime;
            IsConvergent=false;
            while(feCtrlInfo._dt>_MinDt&&!IsConvergent){
                feCtrlInfo._ctan[1]=1.0/feCtrlInfo._dt;
                _CurrentStep+=1;
                _CurrentTime+=feCtrlInfo._dt;
                feCtrlInfo._CurrentStep=_CurrentStep;
                feCtrlInfo._CurrentT=_CurrentTime;
                cout<<"*** Warning: solve failed, dt will be reduced !!!          ***"<<endl;
                if(nonLinearSolver.Solve(mesh,dofHandler,solutionSystem,equationSystem,
                                 elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo)){
                    // if solve success
                    solutionSystem._U=solutionSystem._Unew;
                    solutionSystem.UpdateSolution(_CurrentStep);
                    // update history variables
                    feSystem.FormFE(8,feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,
                            mesh,dofHandler,elmtSystem,mateSystem,
                            solutionSystem._U,solutionSystem._V,
                            solutionSystem._Hist,solutionSystem._HistOld,
                            solutionSystem._ProjValue,
                            equationSystem._AMATRIX,equationSystem._RHS);
                    solutionSystem._HistOld=solutionSystem._Hist;

                    if(_CurrentStep%feCtrlInfo._Interval==0){
                        if(feCtrlInfo._IsDebugOn && !feCtrlInfo._IsDebugDep){
                            nonLinearSolver.PrintIterationInfo();
                        }
                        outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U);
                        printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());
                    }
                    IsConvergent=true;
                    break;
                }
                else{
                    _CurrentStep-=1;
                    _CurrentTime-=feCtrlInfo._dt;
                    feCtrlInfo._CurrentStep=_CurrentStep;
                    feCtrlInfo._CurrentT=_CurrentTime;
                    feCtrlInfo._dt*=_CutbackFactor;
                }
            }
            if(!IsConvergent){
                _TimerEnd=chrono::high_resolution_clock::now();
                _Duration=chrono::duration_cast<std::chrono::microseconds>(_TimerEnd-_TimerStart).count()/1.0e6;
                cout<<"*** Minimal dt is reached, transient analysis(CN) failed!!!***"<<endl;
                cout<<"*** Please check your input file or your elements!!!       ***"<<endl;
                printf("*** Transi analysis finished(fail)     ===>[%12.5f s]***\n",_Duration);
                Msg_AsFem_Exit();
            }
        }
        if(nonLinearSolver.GetCurrentIters()<=_OptIters){
            feCtrlInfo._dt*=_GrowthFactor;
            if(feCtrlInfo._dt>_MaxDt) feCtrlInfo._dt=_MaxDt;
        }
        else{
            feCtrlInfo._dt*=_CutbackFactor;
        }
    }
    _TimerEnd=chrono::high_resolution_clock::now();
    _Duration=chrono::duration_cast<std::chrono::microseconds>(_TimerEnd-_TimerStart).count()/1.0e6;
    cout<<"***--------------------------------------------------------***"<<endl;
    printf("*** Transient analysis finished !      ===>[%12.5f s]***\n",_Duration);
}