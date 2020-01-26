#ifndef ASFEM_NONLINEARSOLVERBLOCKINFO_H
#define ASFEM_NONLINEARSOLVERBLOCKINFO_H

#include <iostream>
#include <iomanip>
#include <string>

#include "NonLinearSolverType.h"

using namespace std;

class NonLinearSolverBlockInfo{
public:
    int _MaxIters=100;
    int _MaxLineSearch=500;
    int _Criterion=0; //0->use either the residual norm or the energy norm
                      //1->use only the residual norm
                      //2->use only the energy norm
                      //3->use the line search ratio norm(equivalent to energy norm) 
    NonLinearSolverType _SolverType=NonLinearSolverType::NEWTONRAPHSON;
    double _Rreltol=1.0e-9,_Rabstol=1.0e-8;
    double _Ereltol=1.0e-21,_Eabstol=1.0e-19;
    double _Stol=1.0e-4;

    void Reset(){
        _MaxIters=100;
        _MaxLineSearch=500;
        _Criterion=0;
        _SolverType=NonLinearSolverType::NEWTONRAPHSON;
        _Rreltol=1.0e-9;_Rabstol=1.0e-8;
        _Ereltol=1.0e-21;_Eabstol=1.0e-19;
        _Stol=1.0e-4;
    }
};

#endif // ASFEM_NONLINEARSOLVERBLOCKINFO_H