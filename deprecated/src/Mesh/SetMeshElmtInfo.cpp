#include "Mesh/Mesh.h"

void Mesh::SetMeshElmtInfo(ElmtSystem &elmtSystem){
    ElmtBlock elmtBlock;
    int iblock,e,i;
    _MeshElmtInfo.resize(GetBulkElmtsNum());

    _BulkElmtDofIndexList.resize(GetBulkElmtsNum(),vector<int>(0));
    
    for(iblock=1;iblock<=elmtSystem.GetElmtBlockNum();++iblock){
        elmtBlock=elmtSystem.GetIthElmtBlock(iblock);
        for(e=1;e<=GetElmtsNumViaPhyName(elmtBlock._ElmtDomainBlockName);++e){
            i=GetElmtIndexViaPhyName(elmtBlock._ElmtDomainBlockName,e);// this is global id(bulk+bc)
            _MeshElmtInfo[i-GetBCElmtsNum()-1].first=elmtBlock._ElmtTypeID;
            _MeshElmtInfo[i-GetBCElmtsNum()-1].second.first=elmtBlock._MateTypeID;
            _MeshElmtInfo[i-GetBCElmtsNum()-1].second.second=iblock;
            _BulkElmtDofIndexList[i-GetBCElmtsNum()-1]=elmtBlock._ElmtDofIndexList;
        }
    }
}