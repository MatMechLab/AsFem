#ifndef ASFEM_NONLINEARSOLVER_H
#define ASFEM_NONLINEARSOLVER_H

#include <iostream>
#include <iomanip>

#include "LinearSolver/LinearSolver.h"
#include "LinearSolver/LinearSolverType.h"
#include "NonLinearSolverType.h"
#include "MessagePrint/MessagePrint.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "BCs/BCSystem.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MaterialSystem/MaterialSystem.h"
#include "FESystem/FESystem.h"

#include "FEProblem/FECtrlInfo.h"

#include "LinearSolver/LinearSolverBlockInfo.h"
#include "NonLinearSolverBlockInfo.h"

class NonLinearSolver
{
public:
    NonLinearSolver(NonLinearSolverType method=NonLinearSolverType::NEWTONRAPHSON);
    void InitNonlinearSolver(NonLinearSolverBlockInfo &nonlinearSolverBlockInfo,
                             LinearSolverBlockInfo &linearSolverBlockInfo);
    void InitLinearSolver(Eigen::SparseMatrix<double,Eigen::RowMajor> &A,Eigen::VectorXd &F,Eigen::VectorXd &x){
        _linearSolver.Solve(A,F,x);
    }

    bool Solve(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);

    // For linear solver
    void SetLinearSolverType(LinearSolverType type) {_linearSolver.SetSolverType(type);}
    void SetLinearSolverMaxIters(const int &maxiters) {_linearSolver.SetMaxIters(maxiters);}
    void SetLinearSolverTolerance(const double &tol) {_linearSolver.SetTolerance(tol);}

    // for Nonlinear method
    void SetNonLinearSolverType(NonLinearSolverType type){_SolverType=type;}
    
    // settings to nonlinear solver system
    void SetDebugMode(bool isdebug=false) {_IsDebug=isdebug;_IsDepDebug=false;}
    void SetDepDebugMode(bool isdebug=false) {_IsDepDebug=isdebug;}
    void SetMaxIters(int maxiters=15) {_MaxIters=maxiters;}
    void SetResidualRelativeError(double rtol=1.0e-8){_RRelTol=rtol;}
    void SetResidualAbsoluteError(double atol=1.0e-8){_RAbsTol=atol;}
    void SetEnergyRelativeError(double rtol=1.0e-18){_ERelTol=rtol;}
    void SetEnergyAbsoluteError(double atol=1.0e-18){_EAbsTol=atol;}
    void SetConvergenceCriterion(int mode=0) {_ConvergenceCriterion=mode;}
    // get information of nonlinear solver
    inline NonLinearSolverType GetNonlinearSolverType() const {return _SolverType;}
    inline int GetCurrentIters() const {return _CurrentIters;}
    inline int GetMaxIters() const {return _MaxIters;}
    inline double GetResidualRelativeError() const {return _RRelTol;}
    inline double GetResidualAbsoluteError() const {return _RAbsTol;}
    inline double GetEnergyRelativeError() const {return _ERelTol;}
    inline double GetEnergyAbsoluteError() const {return _EAbsTol;}
    
    // very important information for adaptive time stepping
    inline double GetResidualInitNorm() const {return _Rnorm0;}
    inline double GetResidualNorm() const {return _Rnorm;}
    inline double GetEnergyInitNorm() const {return _Enorm0;}
    inline double GetEnergyNorm() const {return _Enorm;}

    inline int GetConvergenceCriterion() const {return _ConvergenceCriterion;}

    void PrintNonLinearSolverInfo() const;

    void PrintDepIterationInfo() const;
    void PrintIterationInfo() const;

private:
    bool Linear(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool NewtonRaphson(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool NewtonRaphsonWithLineSearch(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool NewtonRaphsonWithBisectionLineSearch(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool NewtonRaphsonWithRegulaFalsiLineSearch(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool NewtonRaphsonModify(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);
    bool BFGS(Mesh &mesh,DofHandler &dofHandler,
               SolutionSystem &solutionSystem,EquationSystem &equationSystem,
               ElmtSystem &elmtSystem,MaterialSystem &mateSystem,BCSystem &bcSystem,
               FESystem &feSystem,
               FECtrlInfo &feCtrlInfo);

    bool ConvergenceCheck();


    void PrintIterationProcess();


private:
    bool _IsDebug=false,_IsDepDebug=false;
    int _ConvergenceCriterion;
    NonLinearSolverType _SolverType;
    int _MaxIters,_CurrentIters;
    double _RAbsTol,_RRelTol;
    double _EAbsTol,_ERelTol;
    double _Rnorm0,_Rnorm;
    double _Enorm0,_Enorm;
    double _Stol=0.8;

    // for line search
    double _s=1.0,_eta=1.0,_eta0=1.0,_s0=1.0;
    int _MaxLineSearch=100,_nLineSearch=0;
    double _sL,_sU,_etaL,_etaU;

    LinearSolver _linearSolver;

};

#endif // ASFEM_NONLINEARSOLVER_H