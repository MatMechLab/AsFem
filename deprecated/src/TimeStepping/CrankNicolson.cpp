#include "TimeStepping/TimeStepping.h"

void TimeStepping::CrankNicolson(Mesh &mesh,DofHandler &dofHandler,
        SolutionSystem &solutionSystem,EquationSystem &equationSystem,
        ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
        BCSystem &bcSystem,ICSystem &icSystem,
        FESystem &feSystem,NonLinearSolver &nonLinearSolver,
        FECtrlInfo &feCtrlInfo,
        OutputSystem &outputSystem){
    // init the basic settings
    _dt0=feCtrlInfo._dt0;
    _dt=feCtrlInfo._dt0;
    _CurrentTime=0.0;
    _TotalTime=feCtrlInfo._EndTime;

    feCtrlInfo._TimeSteppingMethod=TimeSteppingType::CRANKNICOLSON;
    feCtrlInfo._ISW=6;

    feCtrlInfo._ctan[0]=0.5;
    feCtrlInfo._ctan[1]=1.0/feCtrlInfo._dt0;

    icSystem.ApplyIC(mesh,dofHandler,solutionSystem._U);
    solutionSystem._Unew=solutionSystem._U;


    _TotalStep=(int)(feCtrlInfo._EndTime/feCtrlInfo._dt0);
    _CurrentStep=0;
    outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U);
    printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());

    for(_CurrentTime=0.0;_CurrentTime<=_TotalTime;){
        if(feCtrlInfo._IsDebugOn || feCtrlInfo._IsDebugDep){
                printf("*** Step=%7d, t=%12.6e, dt=%12.6e          ***\n",(int)_CurrentStep,_CurrentTime,feCtrlInfo._dt);
            }
        if(nonLinearSolver.Solve(mesh,dofHandler,solutionSystem,equationSystem,
                                 elmtSystem,mateSystem,bcSystem,feSystem,feCtrlInfo)){
            _CurrentStep+=1;
            _CurrentTime+=feCtrlInfo._dt;
            solutionSystem._U=solutionSystem._Unew;
            solutionSystem.UpdateSolution(_CurrentStep);
            if(_CurrentStep%feCtrlInfo._Interval==0){
                if(feCtrlInfo._IsDebugOn && !feCtrlInfo._IsDebugDep){
                    nonLinearSolver.PrintIterationInfo();
                }
                outputSystem.WriteResultToVTU(_CurrentStep,mesh,dofHandler,solutionSystem._U);
                printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());
            }
        }
        else{
            _TimerEnd=chrono::high_resolution_clock::now();
            _Duration=chrono::duration_cast<std::chrono::microseconds>(_TimerEnd-_TimerStart).count()/1.0e6;
            cout<<"*** Something is wrong, transient analysis(CN) failed!     ***"<<endl;
            cout<<"*** Please check your input file or your elements!!!       ***"<<endl;
            printf("*** Transi analysis finished(fail)     ===>[%12.5f s]***\n",_Duration);
            Msg_AsFem_Exit();
        }
    }
    _TimerEnd=chrono::high_resolution_clock::now();
    _Duration=chrono::duration_cast<std::chrono::microseconds>(_TimerEnd-_TimerStart).count()/1.0e6;
    cout<<"***--------------------------------------------------------***"<<endl;
    printf("*** Transient analysis finished !      ===>[%12.5f s]***\n",_Duration);
}