#include "DofHandler/DofHandler.h"

void DofHandler::ModifyDofActiveMapViaBC(Mesh &mesh,BCSystem &bcSystem){
    // for the nodal dof index with DirichletBC, the value should be zeros
    // otherwise, it is always 1
    int iblock;
    int DofIndex;
    string bcname,bctype;
    BCBlock bcBlock;
    int i,j,e,ee,iInd,nDim;
    int nNodesPerBCElmt=0;
    nDim=mesh.GetDim();
    if(nDim==1){
        nNodesPerBCElmt=1;
    }
    else if(nDim==2){
        nNodesPerBCElmt=mesh.GetNodesNumPerLineElmt();
    }
    else if(nDim==3){
        nNodesPerBCElmt=mesh.GetNodesNumPerSurfaceElmt();
    }
    for(iblock=1;iblock<=bcSystem.GetBCBlockNum();++iblock){
        bcBlock=bcSystem.GetIthBCBlock(iblock);
        DofIndex=GetDofIndexViaName(bcBlock._DofName);
        
        bctype=bcBlock._BCElmtName;
        bcname=bcBlock._BoundaryName;
        
        if(bctype=="dirichlet"){
            for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
                ee=mesh.GetElmtIndexViaPhyName(bcname,e);
                for(i=1;i<=nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=GetIthNodeJthDofIndex(j,DofIndex)-1;// nodal's dof index(global one)
                    // now we need to find the related local one
                    _NodalDofActiveFlag[j-1][iInd]=0.0;
                }
            }
        }
        else {
            continue;
        }
    }
}