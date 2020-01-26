#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem()
{
    _nDofs=0;
    _ZeroCoeffList.clear();
    _CoeffList.clear();
    
    _RhsNorms.clear();
    _dUNorms.clear();
    
    _IsInit=false;
}