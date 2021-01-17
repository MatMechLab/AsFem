//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.18
//+++ Purpose: Create the dofs map for the bulk element
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

void BulkDofHandler::CreateBulkDofsMap(const Mesh &mesh,BCSystem &bcSystem,ElmtSystem &elmtSystem){
    _nMinDim=mesh.GetBulkMeshMinDim();
    _nMaxDim=mesh.GetBulkMeshDim();

    _nNodes=mesh.GetBulkMeshNodesNum();
    _nElmts=mesh.GetBulkMeshElmtsNum();// the total number of element(include both the bulk and bc elmts)
    _nBulkElmts=mesh.GetBulkMeshBulkElmtsNum();

    _nNodesPerBulkElmt=mesh.GetBulkMeshNodesNumPerBulkElmt();
    _nMaxDofsPerElmt=_nNodesPerBulkElmt*_nDofsPerNode;

    _NodalDofFlag.resize(_nNodes,vector<double>(_nDofsPerNode,-1.0));// for dirichlet bc is  : 0
                                                                    // for neumann bc is    : 1
                                                                    // for robin/other bc is: 2
    _ElmtDofFlag.resize(_nBulkElmts,vector<double>(_nMaxDofsPerElmt,-1.0));

    _NodeDofsMap.resize(_nNodes,vector<int>(_nDofsPerNode,-1));
    _ElmtDofsMap.resize(_nBulkElmts,vector<int>(_nMaxDofsPerElmt,-1));

    vector<pair<ElmtType,MateType>> temp;
    temp.clear();
    _ElmtElmtMateTypePairList.resize(_nBulkElmts,temp);
    for(auto &it:_ElmtElmtMateTypePairList){
        it.clear();
    }
    _ElmtLocalDofIndex.resize(_nBulkElmts,vector<vector<int>>(0));
    _ElmtElmtMateIndexList.resize(_nBulkElmts,vector<int>(0));


    _nDofs=_nNodes*_nDofsPerNode;

    _nActiveDofs=0;

    
    int iblock,e,ee,i,ii,j,k,iInd,jInd,ndofs;
    string domainname;
    vector<int> dofindex;

    string str;
    char buff[50];



    //*** firstly, we check whether all the domain has been assigned an [elmt] block
    bool DomainHasElmtBlock=false;
    for(i=1;i<=mesh.GetBulkMeshPhysicalGroupNum();i++){
        DomainHasElmtBlock=false;
        if(mesh.GetBulkMeshIthPhysicalDim(i)==_nMaxDim){// only check the bulk elmts
            DomainHasElmtBlock=false;
            for(j=1;j<=elmtSystem.GetBulkElmtBlockNums();j++){
                cout<<"j="<<j<<": "<<mesh.GetBulkMeshIthPhysicalName(i)<<", "<<elmtSystem.GetIthBulkElmtBlock(j)._DomainName<<endl;
                if(mesh.GetBulkMeshIthPhysicalName(i)==elmtSystem.GetIthBulkElmtBlock(j)._DomainName){
                    DomainHasElmtBlock=true;
                    break;
                }
            }
            if(!DomainHasElmtBlock){
                break;
            }
        }
    }
    if(!DomainHasElmtBlock){
        snprintf(buff,50,"domain %s has no elmt block assigned!",mesh.GetBulkMeshIthPhysicalName(i).c_str());
        str=buff;
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::PrintErrorTxt("please check your input file carefully, all the domain should be assigned at least one [elmt] block!");
        MessagePrinter::AsFem_Exit();
    }

    ElmtType elmttype;
    MateType matetype;
    int mateindex;
    for(iblock=1;iblock<=elmtSystem.GetBulkElmtBlockNums();iblock++){
        domainname=elmtSystem.GetIthBulkElmtBlock(iblock)._DomainName;
        dofindex=elmtSystem.GetIthBulkElmtBlock(iblock)._DofsIDList; // the dof index could be discontinue case, i.e. 1,2,4 !!!
        ndofs=elmtSystem.GetIthBulkElmtBlock(iblock)._nDofs; // the total dofs of current elmt block
        elmttype=elmtSystem.GetIthBulkElmtBlock(iblock)._ElmtType;
        matetype=elmtSystem.GetIthBulkElmtBlock(iblock)._MateType;
        mateindex=elmtSystem.GetIthBulkElmtBlock(iblock)._MateIndex;
        for(auto e:mesh.GetBulkMeshElmtIDsViaPhysicalName(domainname)){
            // now we are in the elmt id vector
            ee=e-(mesh.GetBulkMeshElmtsNum()-mesh.GetBulkMeshBulkElmtsNum());
            _ElmtElmtMateTypePairList[ee-1].push_back(make_pair(elmttype,matetype));
            _ElmtLocalDofIndex[ee-1].push_back(dofindex);
            _ElmtElmtMateIndexList[ee-1].push_back(mateindex);
            for(i=1;i<=mesh.GetBulkMeshIthBulkElmtNodesNum(ee);i++){
                iInd=mesh.GetBulkMeshIthBulkElmtJthNodeID(ee,i);
                for(j=1;j<=ndofs;j++){
                    jInd=dofindex[j-1];
                    _NodeDofsMap[iInd-1][jInd-1]=1;
                    _NodalDofFlag[iInd-1][jInd-1]=1.0;
                }
            }
        }
    }

    // now we can account for the active dofs
    _nActiveDofs=0;
    for(i=1;i<=_nNodes;i++){
        for(j=1;j<=_nDofsPerNode;j++){
            if(_NodalDofFlag[i-1][j-1]>0.0){
                _nActiveDofs+=1;
                _NodeDofsMap[i-1][j-1]=_nActiveDofs;
            }
        }
    }

    // in this case, we reset the bc nodal's dof flag to be zero or other values according to their bc type
    string bcname;
    BCBlock bcBlock;
    for(iblock=1;iblock<=bcSystem.GetBCBlockNums();iblock++){
        bcBlock=bcSystem.GetIthBCBlock(iblock);
        j=bcBlock._DofID;
        if(bcBlock._BCType==BCType::DIRICHLETBC){
            for(auto bc:bcBlock._BoundaryNameList){
                for(auto e:mesh.GetBulkMeshElmtIDsViaPhysicalName(bc)){
                    // now we are in the elmt id vector
                    for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNum(e);i++){
                        iInd=mesh.GetBulkMeshIthElmtJthNodeID(e,i);
                        _NodalDofFlag[iInd-1][j-1]=0.0;
                    }
                }
            }
        }
        else if(bcBlock._BCType==BCType::NODALDIRICHLETBC){
            for(auto bc:bcBlock._BoundaryNameList){
                for(auto i:mesh.GetBulkMeshNodeIDsViaPhysicalName(bc)){
                    _NodalDofFlag[i-1][j-1]=0.0;
                }
            }
        }
        else{
            continue;
        }
    }

    // now we remove all the empty space of some vectors
    _RowNNZ.resize(_nActiveDofs,0);
    _RowMaxNNZ=0;
    for(e=1;e<=_nBulkElmts;e++){
        for(j=1;j<=mesh.GetBulkMeshIthBulkElmtNodesNum(e);j++){
            iInd=mesh.GetBulkMeshIthBulkElmtJthNodeID(e,j);
            for(k=1;k<=_nDofsPerNode;k++){
                ii=(j-1)*_nDofsPerNode+k-1;

                if(_NodalDofFlag[iInd-1][k-1]>=0.0){
                    _ElmtDofsMap[e-1][ii]=_NodeDofsMap[iInd-1][k-1];
                    _RowNNZ[_NodeDofsMap[iInd-1][k-1]-1]+=_nMaxDofsPerElmt;

                    if(_RowNNZ[_NodeDofsMap[iInd-1][k-1]-1]>_RowMaxNNZ){
                        _RowMaxNNZ=_RowNNZ[_NodeDofsMap[iInd-1][k-1]-1];
                    }
                    if(_NodalDofFlag[iInd-1][k-1]>0.0){
                        _ElmtDofFlag[e-1][ii]=1.0;
                    }
                    else if(_NodalDofFlag[iInd-1][k-1]==0.0){
                        _ElmtDofFlag[e-1][ii]=0.0;
                    }
                }
            }
        }
        // remove the zero element space
        _ElmtDofFlag[e-1].erase(
            remove(_ElmtDofFlag[e-1].begin(),_ElmtDofFlag[e-1].end(),-1.0),
            _ElmtDofFlag[e-1].end()
        );
        _ElmtDofFlag[e-1].shrink_to_fit();

        _ElmtDofsMap[e-1].erase(
            remove(_ElmtDofsMap[e-1].begin(),_ElmtDofsMap[e-1].end(),-1),
            _ElmtDofsMap[e-1].end()
        );
        _ElmtDofsMap[e-1].shrink_to_fit();
    }

}