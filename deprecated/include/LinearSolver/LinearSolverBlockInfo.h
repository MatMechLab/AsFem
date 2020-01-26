#ifndef ASFEM_LINEARSOLVERBLOCKINFO_H
#define ASFEM_LINEARSOLVERBLOCKINFO_H

#include <iostream>
#include <iomanip>
#include <string>

#include "LinearSolverType.h"

using namespace std;


class LinearSolverBlockInfo{
public:
    int _MaxIters=2000;// for iterative solver
    double _Tol=1.0e-6;    // for iterative solver
    int _Restart=500;
    string _SolverName="lu";
    LinearSolverType _SolverType=LinearSolverType::SPARSELU;

    void Reset(){
        _MaxIters=2000;// for iterative solver
        _Tol=1.0e-6;    // for iterative solver
        _Restart=500;
        _SolverName="lu";
        _SolverType=LinearSolverType::SPARSELU;
    }
};

#endif // 