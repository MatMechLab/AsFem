#include "FEProblem/FEProblem.h"

void FEProblem::InitFEProblem(){

    _TimerStartOfFEProInit=chrono::high_resolution_clock::now();
    cout<<"*** Start to initialize FE problem ...                     ***"<<endl;


    _TimerStartOfElmtInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize the mesh elmt info ...           ***"<<endl;
    mesh.SetMeshElmtInfo(elmtSystem);
    bcSystem.Init(mesh);
    _TimerEndOfElmtInit=chrono::high_resolution_clock::now();
    _DurationOfElmtInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfElmtInit-_TimerStartOfElmtInit).count()/1.0e6;
    printf("***   mesh elmt info initialized!      ===>[%12.5f s]***\n",_DurationOfElmtInit);
    


    _TimerStartOfDofInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize the dof system ...               ***"<<endl;
    dofHandler.Init(mesh);
    dofHandler.CreateDofMap(mesh,bcSystem);
    // dofHandler.ModifyDofActiveMapViaBC(mesh,bcSystem);
    _TimerEndOfDofInit=chrono::high_resolution_clock::now();
    _DurationOfDofInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfDofInit-_TimerStartOfDofInit).count()/1.0e6;
    printf("***   dof system initialized!          ===>[%12.5f s]***\n",_DurationOfDofInit);
    

    _TimerStartOfSparseInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize the equation system ...          ***"<<endl;
    equationSystem.SetDofsNum(dofHandler.GetActiveDofsNum());
    equationSystem.Init();
    equationSystem.CreateSparsityPatterns(dofHandler);
    _TimerEndOfSparseInit=chrono::high_resolution_clock::now();
    _DurationOfSparseInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfSparseInit-_TimerStartOfSparseInit).count()/1.0e6;
    printf("***   equation system initialized!     ===>[%12.5f s]***\n",_DurationOfSparseInit);

    


    _TimerStartOfFESysInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize FE system ...                    ***"<<endl;
    feSystem.InitFESystem(mesh,dofHandler,qpBlockInfo);
    materialSystem.SetMateValNums(100);
    materialSystem.InitMateValues();//allocate memory for material

    outputSystem.SetInputFileName(inputSystem.GetInputFileName());
    outputSystem.InitOutputStream();

    _TimerEndOfFESysInit=chrono::high_resolution_clock::now();
    _DurationOfFESysInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfFESysInit-_TimerStartOfFESysInit).count()/1.0e6;
    printf("***   FE system initialized!           ===>[%12.5f s]***\n",_DurationOfFESysInit);
    
    _TimerStartOfSolInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize the solution system ...          ***"<<endl;
    solutionSystem.SetDofsNum(dofHandler.GetActiveDofsNum());
    solutionSystem.SetElmtsNum(mesh.GetBulkElmtsNum());
    solutionSystem.SetNodesNum(mesh.GetNodesNum());
    solutionSystem.SetHistsNum(5);// number of hist values on each gauss point
    solutionSystem.SetGPointsNumPerBulkElmt(feSystem.GetBulkElmtGPointNums());
    solutionSystem.Init();
    // must do the following, otherwise the whole system is wrong
    feSystem.SetHistNumPerGPoint(solutionSystem.GetHistNumPerGPoint());
    feSystem.SetProjNumPerNode(solutionSystem.GetProjNumPerNode());
    _TimerEndOfSolInit=chrono::high_resolution_clock::now();
    _DurationOfSolInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfSolInit-_TimerStartOfSolInit).count()/1.0e6;
    printf("***   solution system initialized!     ===>[%12.5f s]***\n",_DurationOfSolInit);

    _TimerStartOfSolverInit=chrono::high_resolution_clock::now();
    cout<<"***  +start to initialize solver system ...                ***"<<endl;
    nonlinearSolver.SetDebugMode(feCtrlInfo._IsDebugOn);
    nonlinearSolver.SetDepDebugMode(feCtrlInfo._IsDebugDep);
    nonlinearSolver.InitNonlinearSolver(nonlinearSolverBlock,linearSolverBlock);
   
    nonlinearSolver.InitLinearSolver(equationSystem._AMATRIX,equationSystem._RHS,equationSystem._dU);
    timeStepping.InitTimeStepping();
    timeStepping.SetAdaptiveFlag(feCtrlInfo._IsAdaptive);
    timeStepping.SetTimeSteppingMethod(feCtrlInfo._TimeSteppingMethod);
    timeStepping.SetMaxDt(feCtrlInfo._DtMax);
    timeStepping.SetMinDt(feCtrlInfo._DtMin);
    timeStepping.SetOptimIters(feCtrlInfo._OptIters);
    _TimerEndOfSolverInit=chrono::high_resolution_clock::now();
    _DurationOfSolverInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfSolverInit-_TimerStartOfSolverInit).count()/1.0e6;
    printf("***   solver system initialized!       ===>[%12.5f s]***\n",_DurationOfSolverInit);



    _TimerEndOfFEProInit=chrono::high_resolution_clock::now();
    _DurationOfFEProInit=chrono::duration_cast<std::chrono::microseconds>(_TimerEndOfFEProInit-_TimerStartOfFEProInit).count()/1.0e6;
    printf("*** FE problem init finished!          ===>[%12.5f s]***\n",_DurationOfFEProInit);


    if(feCtrlInfo._IsPrintMesh){
        mesh.PrintMeshDetailInfo();
    }
    else{
        mesh.PrintMeshInfo();
    }
    dofHandler.PrintDofInfo();
    // dofHandler.PrintDofMap();
    elmtSystem.PrintElmtBlockInfo();
    materialSystem.PrintMateBlockInfo();
    
    
    bcSystem.PrintBCBlockInfo();
    icSystem.PrintICBlockInfo();
    nonlinearSolver.PrintNonLinearSolverInfo();

    JobType=feCtrlInfo._JobType;

}