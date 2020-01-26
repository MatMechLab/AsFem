#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::UpdateSolution(const int &currentstep){
    if(currentstep==1){
        // for the first step
        _Uold.setZero();
        _Uolder.setZero();

        _Vold.setZero();
        _Volder.setZero();

        _Aold.setZero();
        _Aolder.setZero();

        _UNorms.clear();
        _UNorms.push_back(_U.norm());
    }
    else if(currentstep==2){
        _Uold=_U;
        _Uolder.setZero();

        _Vold=_V;
        _Volder.setZero();

        _Aold=_A;
        _Aolder.setZero();

        _UNorms.push_back(_U.norm());
    }
    else{
        _Uolder=_Uold;
        _Uold=_U;

        _Volder=_Vold;
        _Vold=_V;

        _Aolder=_Aold;
        _Aold=_A;

        _UNorms.push_back(_U.norm());
    }
}