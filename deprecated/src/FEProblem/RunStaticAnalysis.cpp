#include "FEProblem/FEProblem.h"

void FEProblem::RunStaticAnalysis(){
    _TimerStartOfStaticJob=chrono::high_resolution_clock::now();
    cout<<"*** Start to run static analysis ...                       ***"<<endl;

    feCtrlInfo._ISW=6;
    feCtrlInfo._CurrentT=0.0;feCtrlInfo._CurrentStep=0;feCtrlInfo._CurrentDt=0.0;
    if(nonlinearSolver.Solve(mesh,dofHandler,solutionSystem,equationSystem,
                             elmtSystem,materialSystem,bcSystem,feSystem,feCtrlInfo)){

        if(feCtrlInfo._IsDebugOn&&!feCtrlInfo._IsDebugDep){
            nonlinearSolver.PrintIterationInfo();
        }
        solutionSystem._U=solutionSystem._Unew;
        _TimerEndOfStaticJob=chrono::high_resolution_clock::now();
        _DurationOfStaticJob=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfStaticJob-_TimerStartOfStaticJob).count()/1.0e6;
        printf("*** Static analysis finished!          ===>[%12.5f s]***\n",_DurationOfStaticJob);
        if(feCtrlInfo._IsProjOn){
            feSystem.FormFE(9,feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,
                            mesh,dofHandler,elmtSystem,materialSystem,
                            solutionSystem._U,solutionSystem._V,
                            solutionSystem._Hist,solutionSystem._HistOld,
                            solutionSystem._ProjValue,
                            equationSystem._AMATRIX,equationSystem._RHS);
            outputSystem.WriteResultToVTU(mesh,dofHandler,solutionSystem._U,
                                          solutionSystem.GetProjNumPerNode(),
                                          solutionSystem.GetProjNameList(),
                                          solutionSystem._ProjValue);
        }
        else{
            outputSystem.WriteResultToVTU(mesh,dofHandler,solutionSystem._U);
        }
        cout<<"***--------------------------------------------------------***"<<endl;
        printf("*** Write results to [%35s] ***\n",outputSystem.GetVTUFileName().c_str());
    }
    else{
        _TimerEndOfStaticJob=chrono::high_resolution_clock::now();
        _DurationOfStaticJob=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfStaticJob-_TimerStartOfStaticJob).count()/1.0e6;
        cout<<"*** Something is wrong, static analysis failed!!!          ***"<<endl;
        cout<<"*** Please check your input file or your elements!!!       ***"<<endl;
        printf("*** Static analysis finished(fail)     ===>[%12.5f s]***\n",_DurationOfStaticJob);
    }
}