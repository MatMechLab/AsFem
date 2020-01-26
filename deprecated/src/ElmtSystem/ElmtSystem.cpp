#include "ElmtSystem/ElmtSystem.h"


ElmtSystem::ElmtSystem()
{
    _nElmtBlock=0;
    _ElmtBlockList.clear();
    _grad_test.setZero();_grad_phi.setZero();
}