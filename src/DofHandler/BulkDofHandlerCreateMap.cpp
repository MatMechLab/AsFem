//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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

void BulkDofHandler::CreateBulkMeshDofsMap(const Mesh &mesh,BCSystem &bcSystem,ElmtSystem &elmtSystem){
    _nMinDim=mesh.GetBulkMeshMinDim();
    _nMaxDim=mesh.GetBulkMeshDim();

    _nNodes=mesh.GetBulkMeshNodesNum();
    _nNodesPerBulkElmt=mesh.GetBulkMeshNodesNumPerBulkElmt();
    _nElmts=mesh.GetBulkMeshElmtsNum();// the total number of element(include both the bulk and bc elmts)
    _nBulkElmts=mesh.GetBulkMeshBulkElmtsNum();

    _nNodesPerBulkElmt=mesh.GetBulkMeshNodesNumPerBulkElmt();
    _nMaxDofsPerElmt=_nNodesPerBulkElmt*_nDofsPerNode;

    _NodalDofFlag.resize(_nNodes,vector<double>(_nDofsPerNode,-1.0));// for dirichlet bc is  : 0
                                                                     // for neumann bc is    : 1
                                                                     // for robin/other bc is: 2
    _BulkElmtDofFlag.resize(_nBulkElmts,vector<double>(_nMaxDofsPerElmt,-1.0));

    _NodeDofsMap.resize(_nNodes,vector<int>(_nDofsPerNode,-1));
    _BulkElmtDofsMap.resize(_nBulkElmts,vector<int>(_nMaxDofsPerElmt,-1));

    vector<pair<ElmtType,MateType>> temp;
    temp.clear();
    _BulkElmtElmtMateTypePairList.resize(_nBulkElmts,temp);
    for(auto &it:_BulkElmtElmtMateTypePairList){
        it.clear();
    }
    _BulkElmtLocalDofIndex.resize(_nBulkElmts,vector<vector<int>>(0));
    _BulkElmtElmtMateIndexList.resize(_nBulkElmts,vector<int>(0));


    _nDofs=_nNodes*_nDofsPerNode;// this is our total dofs

    _nActiveDofs=0;

    
    int iblock,e,ee,i,ii,j,k,iInd,jInd,ndofs;
    string domainname;
    vector<int> dofindex;

    string str;
    char buff[50];



    //*** firstly, we check whether all the domain has been assigned an [elmt] block
    bool DomainHasElmtBlock=false;
    bool UseAllDomain,UseSpeceficDomain;
    for(i=1;i<=mesh.GetBulkMeshPhysicalGroupNum();i++){
        DomainHasElmtBlock=false;UseAllDomain=false;UseSpeceficDomain=false;
        if(mesh.GetBulkMeshIthPhysicalDim(i)==_nMaxDim){// only check the bulk elmts
            for(j=1;j<=elmtSystem.GetBulkElmtBlockNums();j++){
                if(mesh.GetBulkMeshIthPhysicalName(i)==elmtSystem.GetIthBulkElmtBlock(j)._DomainName){
                    if(elmtSystem.GetIthBulkElmtBlock(j)._DomainName=="alldomain"){
                        UseAllDomain=true;UseSpeceficDomain=false;
                    }
                    else{
                        UseAllDomain=false;UseSpeceficDomain=true;
                    }
                    break;
                }
            }
        }
        if(UseAllDomain||UseSpeceficDomain){
            DomainHasElmtBlock=true;
            break;
        }
    }
    if(!DomainHasElmtBlock){
        snprintf(buff,50,"domain %s has no elmt block assigned!",mesh.GetBulkMeshIthPhysicalName(i).c_str());
        str=buff;
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::PrintErrorTxt("please check your input file carefully, all the domain should be assigned at least one [elmt] block!");
        MessagePrinter::AsFem_Exit();
    }

    //*** secondly, we check whether all the DoFs have been assigned by an [elmt] block
    bool DofHasElmtBlock=false;
    string dofname;
    for(i=1;i<=GetDofsNumPerNode();i++){
        DofHasElmtBlock=false;
        dofname=GetIthDofName(i);
        for(j=1;j<=elmtSystem.GetBulkElmtBlockNums();j++){
            for(const auto &name:elmtSystem.GetIthBulkElmtBlock(j)._DofsNameList){
                if(name==dofname){
                    DofHasElmtBlock=true;
                    break;
                }
            }
            if(DofHasElmtBlock){
                break;
            }
        }
        if(!DofHasElmtBlock){
            break;
        }
    }
    if(!DofHasElmtBlock){
        str="DoF->"+dofname+" has not been assigned to any [elmts] block! Please check your input file carefully";
        MessagePrinter::PrintErrorTxt(str);
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
            _BulkElmtElmtMateTypePairList[ee-1].push_back(make_pair(elmttype,matetype));
            _BulkElmtLocalDofIndex[ee-1].push_back(dofindex);
            _BulkElmtElmtMateIndexList[ee-1].push_back(mateindex);
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
        if(bcBlock._BCType==BCType::DIRICHLETBC||
           bcBlock._BCType==BCType::USER1DIRICHLETBC||
           bcBlock._BCType==BCType::USER2DIRICHLETBC||
           bcBlock._BCType==BCType::USER3DIRICHLETBC||
           bcBlock._BCType==BCType::USER4DIRICHLETBC||
           bcBlock._BCType==BCType::USER5DIRICHLETBC){
            for(auto bc:bcBlock._BoundaryNameList){
                for(auto e:mesh.GetBulkMeshElmtIDsViaPhysicalName(bc)){
                    // now we are in the elmt id vector
                    for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNum(e);i++){
                        iInd=mesh.GetBulkMeshIthElmtJthNodeID(e,i);
                        for(const auto &dofid:bcBlock._DofIDs){
                            _NodalDofFlag[iInd-1][dofid-1]=0.0;
                        }
                    }
                }
            }
        }
        else if(bcBlock._BCType==BCType::NODALDIRICHLETBC){
            for(auto bc:bcBlock._BoundaryNameList){
                for(auto i:mesh.GetBulkMeshNodeIDsViaPhysicalName(bc)){
                    for(const auto &dofid:bcBlock._DofIDs){
                        _NodalDofFlag[i-1][dofid-1]=0.0;
                    }
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
                    _BulkElmtDofsMap[e-1][ii]=_NodeDofsMap[iInd-1][k-1];
                    _RowNNZ[_NodeDofsMap[iInd-1][k-1]-1]+=_nMaxDofsPerElmt;

                    if(_RowNNZ[_NodeDofsMap[iInd-1][k-1]-1]>_RowMaxNNZ){
                        _RowMaxNNZ=_RowNNZ[_NodeDofsMap[iInd-1][k-1]-1];
                    }
                    if(_NodalDofFlag[iInd-1][k-1]>0.0){
                        _BulkElmtDofFlag[e-1][ii]=1.0;
                    }
                    else if(_NodalDofFlag[iInd-1][k-1]==0.0){
                        _BulkElmtDofFlag[e-1][ii]=0.0;
                    }
                }
            }
        }
        // remove the zero element space
        _BulkElmtDofFlag[e-1].erase(
            remove(_BulkElmtDofFlag[e-1].begin(),_BulkElmtDofFlag[e-1].end(),-1.0),
            _BulkElmtDofFlag[e-1].end()
        );
        _BulkElmtDofFlag[e-1].shrink_to_fit();

        _BulkElmtDofsMap[e-1].erase(
            remove(_BulkElmtDofsMap[e-1].begin(),_BulkElmtDofsMap[e-1].end(),-1),
            _BulkElmtDofsMap[e-1].end()
        );
        _BulkElmtDofsMap[e-1].shrink_to_fit();
    }

}
