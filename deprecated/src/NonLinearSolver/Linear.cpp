#include "NonLinearSolver/NonLinearSolver.h"

bool NonLinearSolver::Linear(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo){
    _CurrentIters=0;
    // feSystem.PresetDirichletBC(mesh,dofHandler,bcSystem,solutionSystem._Unew);
    bcSystem.ApplyPreSetDirichletBC(feCtrlInfo._CurrentT,mesh,dofHandler,solutionSystem._Unew);
    solutionSystem._V.setZero();
   
    feSystem.FormFE(feCtrlInfo._ISW,feCtrlInfo._CurrentT,feCtrlInfo._CurrentDt,feCtrlInfo._ctan,mesh,dofHandler,elmtSystem,mateSystem,
    solutionSystem._Unew,solutionSystem._V,solutionSystem._Hist,solutionSystem._HistOld,
    solutionSystem._ProjValue,
    equationSystem._AMATRIX,equationSystem._RHS);
    
    bcSystem.SetMaxKMatrixValue(feSystem.GetMaxAMatrixValue()*1.0e3);

    bcSystem.ApplyBC(mesh,dofHandler,feSystem.GetFEPtr(),
                    feCtrlInfo._CurrentT,feCtrlInfo._ctan,
                    equationSystem._AMATRIX,equationSystem._RHS,
                    solutionSystem._Unew);
    if(_linearSolver.Solve(equationSystem._AMATRIX,equationSystem._RHS,equationSystem._dU)){
        _Rnorm=equationSystem._RHS.norm();
        _Enorm=equationSystem._dU.norm()*_Rnorm;
        solutionSystem._Unew=-equationSystem._dU;
        return true;
    }
    else{
        cout<<"*** Error: solve Ax=F failed in linear/NonlinearSolver!!!"<<endl;
        return false;
    }
    return true;
}