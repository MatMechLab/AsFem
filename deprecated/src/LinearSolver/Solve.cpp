#include "LinearSolver/LinearSolver.h"

bool LinearSolver::Solve(Eigen::SparseMatrix<double,Eigen::RowMajor> &A,
                         Eigen::VectorXd &F,
                         Eigen::VectorXd &x)
{
    if(!_IsInit){
        switch(GetSolverType()){
            case LinearSolverType::SPARSELU:
                //_SparseLUSolver.setPivotThreshold(0.8);
                _SparseLUSolver.analyzePattern(A);
                _IsInit=true;
                /*
                //here we can't get info, it will trigger the assertion errors
                if(_SparseLUSolver.info()!=Eigen::Success){
                    Msg_LinearSolver_AnalyzeFailed();
                    return false;
                }*/
                break;
            case LinearSolverType::SPARSEQR:
                _SparseQRSolver.analyzePattern(A);
                _IsInit=true;
                // if(_SparseQRSolver.info()!=Eigen::Success){
                //     Msg_LinearSolver_AnalyzeFailed();
                //     return false;
                // }
                break;
            case LinearSolverType::BICG:
                _BICGSolver.setTolerance(_Tolerance);
                _BICGSolver.setMaxIterations(_MaxIters);
                _BICGSolver.analyzePattern(A);
                _IsInit=true;
                // if(_BICGSolver.info()!=Eigen::Success){
                //     Msg_LinearSolver_AnalyzeFailed();
                //     return false;
                // }
                break;
            case LinearSolverType::CG:
                _CGSolver.setMaxIterations(_MaxIters);
                _CGSolver.setTolerance(_Tolerance);
                _CGSolver.analyzePattern(A);
                // if(_CGSolver.info()!=Eigen::Success){
                //     Msg_LinearSolver_AnalyzeFailed();
                //     return false;
                // }
                _IsInit=true;
                break;
            case LinearSolverType::GMRES:
                _GMRESSolver.setTolerance(_Tolerance);
                _GMRESSolver.setMaxIterations(_MaxIters);
                _GMRESSolver.set_restart(_Restarts);
                _GMRESSolver.analyzePattern(A);
                _IsInit=true;
                break;
            case LinearSolverType::PARDISOLU:
                // for pardiso, it seems we don't need to do anythings
                _IsInit=true;
                break;
            case LinearSolverType::UMFPACK:
                _IsInit=true;
                break;
            default:
                assert("*** Error: Unsupported linear solver type in Solver.cpp!!!");
                break;
        }
    }
    else if(_IsInit){
        switch(GetSolverType()){
            case LinearSolverType::SPARSELU:
                _SparseLUSolver.factorize(A);
                if(_SparseLUSolver.info()!=Eigen::Success){
                    Msg_LinearSolver_FactorizeFailed();
                    return false;
                }
                x=_SparseLUSolver.solve(F);
                if(_SparseLUSolver.info()==Eigen::Success){
                    return true;
                }
                else
                {
                    Msg_LinearSolver_SolveFailed();
                    return false;
                }
                break;
            case LinearSolverType::SPARSEQR:
                _SparseQRSolver.factorize(A);
                if(_SparseQRSolver.info()!=Eigen::Success){
                    Msg_LinearSolver_FactorizeFailed();
                    return false;
                }
                x=_SparseQRSolver.solve(F);
                if(_SparseQRSolver.info()==Eigen::Success){
                    return true;
                }
                else
                {
                    Msg_LinearSolver_SolveFailed();
                    return false;
                }
                break;
            case LinearSolverType::CG:
                _CGSolver.factorize(A);
                if(_CGSolver.info()!=Eigen::Success){
                    Msg_LinearSolver_FactorizeFailed();
                    return false;
                }
                x=_CGSolver.solve(F);
                if(_CGSolver.info()==Eigen::Success){
                    return true;
                }
                else
                {
                    Msg_LinearSolver_SolveFailed();
                    return false;
                }
                break;
            case LinearSolverType::BICG:
                 _BICGSolver.factorize(A);
                if(_BICGSolver.info()!=Eigen::Success){
                    Msg_LinearSolver_FactorizeFailed();
                    return false;
                }
                x=_BICGSolver.solve(F);
                if(_BICGSolver.info()==Eigen::Success){
                    return true;
                }
                else
                {
                    Msg_LinearSolver_SolveFailed();
                    return false;
                }
                break;
            case LinearSolverType::GMRES:
                 _GMRESSolver.factorize(A);
                x=_GMRESSolver.solve(F);
                if(_GMRESSolver.info()==Eigen::Success){
                    return true;
                }
                else
                {
                    Msg_LinearSolver_SolveFailed();
                    return false;
                }
                break;
            case LinearSolverType::PARDISOLU:
                #ifdef USE_PARDISO
                _PardisoLUSolver.compute(A);
                x=_PardisoLUSolver.solve(F);
                return true;
                #else
                assert("*** Error: in order to use Pardiso, you must enable it in your cmake file!!!");
                #endif
                break;
            case LinearSolverType::UMFPACK:
                #ifdef USE_SUITESPARSE
                //double *null = (double *) NULL;
                void *Symbolic, *Numeric;
                (void) umfpack_di_symbolic (A.rows(),A.rows(),A.outerIndexPtr(),A.innerIndexPtr(),A.valuePtr(),&Symbolic,nullptr,nullptr);
                (void) umfpack_di_numeric (A.outerIndexPtr(),A.innerIndexPtr(),A.valuePtr(),Symbolic,&Numeric,nullptr,nullptr);
                umfpack_di_free_symbolic (&Symbolic);
                (void)umfpack_di_solve (UMFPACK_A,A.outerIndexPtr(),A.innerIndexPtr(),A.valuePtr(),x.data(),F.data(),Numeric,nullptr,nullptr) ;
                umfpack_di_free_numeric(&Numeric);
                #endif
                break;
            default:
                assert("*** Error: Unsupported linear solver type in Solver.cpp!!!");
                break;
        }
    }
    return true;
}