//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2019
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

void Mesh::SetElmtInfo(ElmtSystem &elmtSystem){
    ElmtBlock elmtBlock;
    int iblock,e,i;

    _MeshElmtInfoList.resize(GetBulkElmtsNum());

    _BulkElmtDofIndexList.resize(GetBulkElmtsNum(),vector<int>(0));
    
    for(iblock=1;iblock<=elmtSystem.GetElmtBlocksNum();++iblock){
        elmtBlock=elmtSystem.GetIthElmtBlock(iblock);
        for(e=1;e<=GetElmtsNumViaPhyName(elmtBlock._DomainName);++e){
            i=GetIthElmtIndexViaPhyName(elmtBlock._DomainName,e);// this is global id(bulk+bc)
            // _MeshElmtInfoList[i-GetBCElmtsNum()-1]=make_pair(elmtBlock._ElmtType,elmtBlock._MateType);
            _MeshElmtInfoList[i-GetBCElmtsNum()-1].first=elmtBlock._ElmtType;
            _MeshElmtInfoList[i-GetBCElmtsNum()-1].second.first=elmtBlock._MateType;
            // _MeshElmtInfoList[i-GetBCElmtsNum()-1].second.second=iblock;
            _MeshElmtInfoList[i-GetBCElmtsNum()-1].second.second=elmtBlock._MateBlockIndex;
            _BulkElmtDofIndexList[i-GetBCElmtsNum()-1]=elmtBlock._DofIndexList;
        }
    }
}