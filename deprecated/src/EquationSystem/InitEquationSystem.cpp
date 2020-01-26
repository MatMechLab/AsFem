#include "EquationSystem/EquationSystem.h"

void EquationSystem::Init()
{
    _ZeroCoeffList.clear();
    _CoeffList.clear();

    _RHS=Eigen::VectorXd(_nDofs);_RHS.setZero();
    _dU=Eigen::VectorXd(_nDofs);_dU.setZero();

    _AMATRIX=Eigen::SparseMatrix<double,Eigen::RowMajor>(_nDofs,_nDofs);

    _IsInit=true;
}