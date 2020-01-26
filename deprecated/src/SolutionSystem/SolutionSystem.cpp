#include "SolutionSystem/SolutionSystem.h"

SolutionSystem::SolutionSystem(){
    _nDofs=0;_nNodes=0;_nBulkElmts=0;
    _nDofsPerNode=0;
    _nProjsPerNode=0;
    _nHistValuePerGPoint=0;
    _nGPointPerBulkElmt=0;
    _HasProjNameList=false;
}