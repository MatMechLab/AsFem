#include "DofHandler/DofHandler.h"

void DofHandler::Init(Mesh &mesh)
{
    _nDim=mesh.GetDim();
    _nNodes=mesh.GetNodesNum();
    _nElmts=mesh.GetElmtsNum();
    _nBulkElmts=mesh.GetBulkElmtsNum();
    _nNodesPerBulkElmt=mesh.GetNodesNumPerBulkElmt();
    _nMaxDofsPerElmt=_nDofsPerNode*_nNodesPerBulkElmt;
    
    _NodalDofFlag.resize(_nNodes,vector<int>(_nDofsPerNode,0));
    _ElmtalDofIndex.resize(_nBulkElmts,vector<int>(_nDofsPerNode*_nNodesPerBulkElmt,0));
    _nDofs=_nDofsPerNode*_nNodes;
    _nActiveDofs=0;
}