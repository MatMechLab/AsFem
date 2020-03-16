//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
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


    // int size,rank;
    // int nelmts,estart,eend,esize;
    // MPI_Comm_size(PETSC_COMM_WORLD,&size);
    // MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    // for(iblock=1;iblock<=elmtSystem.GetElmtBlocksNum();++iblock){
    //     elmtBlock=elmtSystem.GetIthElmtBlock(iblock);

    //     nelmts=GetElmtsNumViaPhyName(elmtBlock._DomainName);
    //     esize=nelmts/size;
    //     estart=rank*esize;
    //     eend=estart+esize;
    //     if(rank==size-1) eend=nelmts;

    //     // for(e=1;e<=GetElmtsNumViaPhyName(elmtBlock._DomainName);++e)
    //     for(ee=estart;ee<eend;ee++){
    //         e=ee+1;
    //         i=GetIthElmtIndexViaPhyName(elmtBlock._DomainName,e);// this is global id(bulk+bc)
    //         _MeshElmtInfoList[i-GetBCElmtsNum()-1].first=elmtBlock._ElmtType;
    //         _MeshElmtInfoList[i-GetBCElmtsNum()-1].second.first=elmtBlock._MateType;
    //         _MeshElmtInfoList[i-GetBCElmtsNum()-1].second.second=elmtBlock._MateBlockIndex;
    //         _BulkElmtDofIndexList[i-GetBCElmtsNum()-1]=elmtBlock._DofIndexList;

    //         if(rank!=0){
    //             MPI_Send(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].first,1,MPI_INT,0,11,PETSC_COMM_WORLD);
    //             MPI_Send(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].second.first,1,MPI_INT,0,12,PETSC_COMM_WORLD);
    //             MPI_Send(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].second.second,1,MPI_INT,0,13,PETSC_COMM_WORLD);
    //             MPI_Send(&_BulkElmtDofIndexList[i-GetBCElmtsNum()-1][0],elmtBlock._DofIndexList.size(),MPI_INT,0,2,PETSC_COMM_WORLD);
    //         }
    //     }
    // }
    // // now we receive all the messages from different processors
    // if(rank==0){
    //     for(iblock=1;iblock<=elmtSystem.GetElmtBlocksNum();++iblock){
    //         elmtBlock=elmtSystem.GetIthElmtBlock(iblock);
    //         nelmts=GetElmtsNumViaPhyName(elmtBlock._DomainName);
    //         esize=nelmts/size;
    //         for(int cpu=1;cpu<size;cpu++){
    //             estart=cpu*esize;
    //             eend=estart+esize;
    //             if(cpu==size-1) eend=nelmts;
    //             for(ee=estart;ee<eend;ee++){
    //                 e=ee+1;
    //                 i=GetIthElmtIndexViaPhyName(elmtBlock._DomainName,e);// this is global id(bulk+bc)
    //                 MPI_Recv(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].first,1,MPI_INT,cpu,11,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
    //                 MPI_Recv(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].second.first,1,MPI_INT,cpu,12,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
    //                 MPI_Recv(&_MeshElmtInfoList[i-GetBCElmtsNum()-1].second.second,1,MPI_INT,cpu,13,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
    //                 MPI_Recv(&_BulkElmtDofIndexList[i-GetBCElmtsNum()-1][0],elmtBlock._DofIndexList.size(),MPI_INT,cpu,2,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
    //             }
    //         }
    //     }
    // }
    
}