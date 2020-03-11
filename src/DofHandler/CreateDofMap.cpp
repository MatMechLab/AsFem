//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "DofHandler/DofHandler.h"

void DofHandler::CreateDofMap(Mesh &mesh,BCSystem &bcSystem){
    _nMaxDim=mesh.GetDim();
    _nMinDim=mesh.GetMinDim();

    _nNodes=mesh.GetNodesNum();
    _nElmts=mesh.GetElmtsNum();
    _nBulkElmts=mesh.GetBulkElmtsNum();

    _nNodesPerBulkElmt=mesh.GetNodesNumPerBulkElmt();
    _nMaxDofsPerElmt=_nMaxDofsPerNode*_nNodesPerBulkElmt;

    _nDofs=_nMaxDofsPerNode*_nNodes;

    int iblock;
    int DofIndex;
    string bcname;
    BCType bctype;
    BCBlock bcBlock;
    int i,ii,j,k,e,ee,iInd,nDim;
    int nNodesPerBCElmt=0;


    _NodalDofFlag.resize(_nNodes,vector<PetscInt>(_nMaxDofsPerNode,0));
    _NodalDofActiveFlag.resize(_nNodes,vector<PetscReal>(_nMaxDofsPerNode,0.0));
    _ElmtalDofIndex.resize(_nBulkElmts,vector<PetscInt>(_nMaxDofsPerNode*_nNodesPerBulkElmt,0));
    _ElmtalDofActiveFlag.resize(_nBulkElmts,vector<PetscReal>(_nMaxDofsPerNode*_nNodesPerBulkElmt,-1.0));
    


    _nActiveDofs=0;
    


    vector<PetscInt> elConn,elDofsIndex;
    elConn.resize(mesh.GetNodesNumPerBulkElmt(),0);

    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        mesh.GetIthBulkElmtConn(e,elConn);
        elDofsIndex=mesh.GetIthBulkElmtDofsIndex(e);
        for(i=1;i<=mesh.GetIthBulkElmtNodesNum(e);++i){
            j=elConn[i-1];
            for(k=1;k<=(PetscInt)elDofsIndex.size();++k){
                _NodalDofFlag[j-1][elDofsIndex[k-1]-1]=1;
                _NodalDofActiveFlag[j-1][elDofsIndex[k-1]-1]=1.0;
            }
        }
    }
    

    // now we start to count for the active dofs
    _nActiveDofs=0;
    for(i=0;i<_nNodes;++i){
        for(j=0;j<_nMaxDofsPerNode;++j){
            if(_NodalDofFlag[i][j]!=0){
                _nActiveDofs+=1;
                _NodalDofFlag[i][j]=_nActiveDofs;
                _NodalDofActiveFlag[i][j]=1.0;
            }
            else{
                continue;
            }
        }
    }
    // modify the nodal dof active flags according to related boundary condition
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
    for(iblock=1;iblock<=bcSystem.GetBCBlocksNum();++iblock){
        bcBlock=bcSystem.GetIthBCBlock(iblock);
        bctype=bcSystem.GetIthBCBlock(iblock)._BCType;
        DofIndex=GetDofIndexViaName(bcBlock._DofName);
        //cout<<"bcname="<<bcname<<", dofindex="<<DofIndex<<endl; 
        if(bctype==BCType::DirichletBC){
            for(int ibc=1;ibc<(int)bcSystem.GetIthBCBlock(iblock)._BoundaryNameList.size();++ibc){
                bcname=bcSystem.GetIthBCBlock(iblock)._BoundaryNameList[ibc-1];
                for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
                    ee=mesh.GetIthElmtIndexViaPhyName(bcname,e);
                    for(i=1;i<=nNodesPerBCElmt;++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
                        // iInd=GetIthNodeJthDofIndex(j,DofIndex)-1;// nodal's dof index(global one)
                        // now we need to find the related local one
                        _NodalDofActiveFlag[j-1][DofIndex-1]=0.0;
                        // cout<<"j="<<j<<", iInd="<<iInd<<endl;
                    }
                }
            }
        }
        else {
            continue;
        }
    }

   
    _RowNNZ.resize(_nActiveDofs,0);
    _RowMaxNNZ=-1;
    for(e=0;e<_nBulkElmts;++e){
        for(j=1;j<=mesh.GetIthBulkElmtNodesNum(e+1);++j){
            iInd=mesh.GetIthBulkElmtJthConn(e+1,j);
            for(k=1;k<=GetDofsNumPerNode();++k){
                ii=(j-1)*GetDofsNumPerNode()+k-1;
                
                if(_NodalDofFlag[iInd-1][k-1]!=0){
                    _ElmtalDofIndex[e][ii]=_NodalDofFlag[iInd-1][k-1];
                    _RowNNZ[_NodalDofFlag[iInd-1][k-1]-1]+=_nMaxDofsPerElmt;

                    if(_RowNNZ[_NodalDofFlag[iInd-1][k-1]-1]>_RowMaxNNZ){
                        _RowMaxNNZ=_RowNNZ[_NodalDofFlag[iInd-1][k-1]-1];
                    }
                }
                
                if(_NodalDofActiveFlag[iInd-1][k-1]>0.0){
                    
                    _ElmtalDofActiveFlag[e][ii]=1.0;
                }
                else{
                    _ElmtalDofActiveFlag[e][ii]=0.0;
                }
            }
        }
        _ElmtalDofIndex[e].erase(
                remove(_ElmtalDofIndex[e].begin(),_ElmtalDofIndex[e].end(),0),
                _ElmtalDofIndex[e].end());
        _ElmtalDofIndex[e].shrink_to_fit();

        _ElmtalDofActiveFlag[e].erase(
                remove(_ElmtalDofActiveFlag[e].begin(),_ElmtalDofActiveFlag[e].end(),-1.0),
                _ElmtalDofActiveFlag[e].end());
        _ElmtalDofActiveFlag[e].shrink_to_fit();
    }
}