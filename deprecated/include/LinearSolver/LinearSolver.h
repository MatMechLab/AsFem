#ifndef ASFEM_LINEARSOLVER_H
#define ASFEM_LINEARSOLVER_H

#include <iostream>
#include <iomanip>
#include <string>

#include "Eigen/Sparse"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include "Eigen/IterativeLinearSolvers"
#include "unsupported/Eigen/IterativeSolvers"

#ifdef USE_PARDISO
#define EIGEN_USE_MKL_ALL
#include "Eigen/PardisoSupport"
#endif

#ifdef USE_SUITESPARSE
#include "umfpack.h"
#endif




#include "LinearSolverType.h"
#include "LinearSolverBlockInfo.h"
#include "MessagePrint/MessagePrint.h"

using namespace std;


class LinearSolver
{
public:
    LinearSolver();

    void InitLinearSolver(LinearSolverBlockInfo &linearSolverBlockInfo);
    void SetSolverType(LinearSolverType type){_SolverType=type;}
    void SetMaxIters(const int &maxiters) {_MaxIters=maxiters;}
    void SetTolerance(const double &tol) {_Tolerance=tol;}

    int GetIters() const;
    double GetTolerance() const;
    LinearSolverType GetSolverType() const {return _SolverType;}

    bool Solve(Eigen::SparseMatrix<double,Eigen::RowMajor> &A,
               Eigen::VectorXd &F,Eigen::VectorXd &x);

    void PrintSolverInfo() const;

private:
    LinearSolverType _SolverType;
    int _MaxIters,_Restarts;
    double _Tolerance;

private:
    bool _IsInit=false;
    // define all the solvers
    //Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::COLAMDOrdering<int>> _SparseLUSolver;
    //Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>> _SparseLUSolver;
    Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>> _SparseLUSolver;
    Eigen::SparseQR<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::COLAMDOrdering<int>> _SparseQRSolver;
    // for iterative solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>,
                             Eigen::Lower|Eigen::Upper> _CGSolver;
    
    // it seems bicg can't be paralled when LUT preconditioner is used!
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,
    //                 Eigen::IncompleteLUT<double,int>> _BICGSolver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>> _BICGSolver;

    Eigen::GMRES<Eigen::SparseMatrix<double,Eigen::RowMajor>> _GMRESSolver;// don't use iLUT preconditioner
    
    // when iLUT preconditionier is used, the openmp has not effect at all!!!
    // Eigen::GMRES<Eigen::SparseMatrix<double,Eigen::RowMajor>,
    //              Eigen::IncompleteLUT<double,int>> _GMRESSolver;

    // // for some solver tags from Viennacl
    // viennacl::linalg::gmres_tag _gmres_tag{1.0e-8,500,20};
    // // viennacl::linalg::gmres_tag _gmres_tag;
    // viennacl::compressed_matrix<double> _vl_spmat;
    // viennacl::vector<double> _vl_rhs,_vl_x;
    // viennacl::linalg::chow_patel_tag _chow_patel_ilu_config;
    // use solver from AMGCL


    // for pardiso from intel!
    // in order to use it, you must install intel's icc compiler!!!
    #ifdef USE_PARDISO
    Eigen::PardisoLU<Eigen::SparseMatrix<double>> _PardisoLUSolver;
    #endif

};

#endif // ASFEM_LINEARSOLVER_H