#include "BCs/BCSystem.h"


void BCSystem::Init(Mesh &mesh){
    _nDim=mesh.GetDim();
    if(_nDim==1){
        _nNodesPerBCElmt=1;
    }
    else if(_nDim==2){
        _nNodesPerBCElmt=mesh.GetNodesNumPerLineElmt();
    }
    else if(_nDim==3){
        _nNodesPerBCElmt=mesh.GetNodesNumPerSurfaceElmt();
    }
    _Area=0.0;
    _ElArea=0.0;
    _bcvalue=0.0;
    _DofIndex=1;

     _elNodes.InitNodes(_nNodesPerBCElmt);
}