#ifndef ASFEM_EQUATIONSYSTEM_H
#define ASFEM_EQUATIONSYSTEM_H

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

#include "MessagePrint/MessagePrint.h"

#include "DofHandler/DofHandler.h"

#include "Eigen/Eigen"
#include "Eigen/Core"
#include "Eigen/Sparse"

using namespace std;

typedef Eigen::Triplet<double> T;

class EquationSystem
{
public:
    EquationSystem();

    void Init();
    void CreateSparsityPatterns(DofHandler &dofHandler);

    void SetDofsNum(const int &ndofs) {_nDofs=ndofs;}
    void AddRhsNorm(const double &normval) {_RhsNorms.push_back(normval);}
    void AddDuNrom(const double &normval) {_dUNorms.push_back(normval);}

    inline void FirstZeroAMATRIX() {_AMATRIX.setFromTriplets(_ZeroCoeffList.begin(),_ZeroCoeffList.end());} 
    inline void ZeroAMATRIX(){_AMATRIX.setFromTriplets(_ZeroCoeffList.begin(),_ZeroCoeffList.end());}

    void SetNNZNum(int nnznum) {_NNZNum=nnznum;}

public:
    vector<T> _ZeroCoeffList;
    vector<T> _CoeffList;


public:
    Eigen::VectorXd _RHS;
    Eigen::VectorXd _dU;
    Eigen::SparseMatrix<double,Eigen::RowMajor> _AMATRIX;

private:
    bool _IsInit=false;
    int _nDofs;
    int _NNZNum;
    vector<double> _RhsNorms;// associated with step
    vector<double> _dUNorms;
};

#endif // ASFEM_EQUATIONSYSTEM_H