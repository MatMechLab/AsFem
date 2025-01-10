//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.12
//+++ Purpose: the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"
#include "MPIUtils/MPIDataBus.h"

BulkElmtSystem::BulkElmtSystem(){
    m_ElmtBlockNum=0;
    m_ElmtBlockList.clear();
    m_GlobalBulkElmtsNum=0;
    m_LocalBulkElmtsNum=0;
    m_LocalElemental_ElmtBlockID.clear();
}

void BulkElmtSystem::addElmtBlock2List(ElmtBlock t_elmtblock){
    if(m_ElmtBlockList.size()<1){
        m_ElmtBlockList.push_back(t_elmtblock);
        m_ElmtBlockNum=1;
    }
    else{
        bool IsExist=false;
        for(int i=0;i<static_cast<int>(m_ElmtBlockList.size());i++){
            if(t_elmtblock.m_ElmtBlockName==m_ElmtBlockList[i].m_ElmtBlockName){
                IsExist=true;
                break;
            }
        }
        if(!IsExist){
            m_ElmtBlockList.push_back(t_elmtblock);
            m_ElmtBlockNum+=1;
        }
        else{
            MessagePrinter::printErrorTxt(t_elmtblock.m_ElmtBlockName+" is already exist, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}

void BulkElmtSystem::init(const FECell &t_fecell){
    m_GlobalBulkElmtsNum=t_fecell.getFECellBulkElmtsNum();
    m_LocalBulkElmtsNum=t_fecell.getLocalFECellBulkElmtsNum();

    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        // assign the global elmtblock id vector
        m_GlobalElemental_ElmtBlockID.resize(m_GlobalBulkElmtsNum,vector<int>(0));
        // loop over all the element block list
        for(int i=1;i<=m_ElmtBlockNum;i++){
            for(const auto &phyname:m_ElmtBlockList[i-1].m_DomainNameList){
                for(const auto &elmtid:t_fecell.getFECellBulkElmtIDSetViaPhyName(phyname)){
                    m_GlobalElemental_ElmtBlockID[elmtid-1].push_back(i);
                }
            }
        }
        // send the elmtblockid to each rank based on their partition info
        int rankid;
        m_LocalElemental_ElmtBlockID.clear();
        for(int e=1;e<=t_fecell.getFECellBulkElmtsNum();e++){
            rankid=t_fecell.getFECellIthBulkElmtRankID(e);
            if(rankid==0){
                m_LocalElemental_ElmtBlockID.push_back(m_GlobalElemental_ElmtBlockID[e-1]);
            }
            else{
                MPIDataBus::sentIntegerVectorToOthers(m_GlobalElemental_ElmtBlockID[e-1],rankid*1000+1,rankid);
            }
        }
    }// end-of-master-rank
    else{
        vector<int> idvec;
        for(int e=0;e<t_fecell.getLocalFECellBulkElmtsNum();e++){
            MPIDataBus::receiveIntegerVectorFromMaster(idvec,rank*1000+1);
            m_LocalElemental_ElmtBlockID.push_back(idvec);
        }
    }// end-of-other-ranks
    MPI_Barrier(MPI_COMM_WORLD);
}

void BulkElmtSystem::printBulkElmtSystemInfo()const{
    MessagePrinter::printNormalTxt("Bulk element system information summary");
    MessagePrinter::printNormalTxt(" number of bulk element blocks = "+to_string(getBulkElmtBlocksNum()));
    for(int i=0;i<m_ElmtBlockNum;i++){
        m_ElmtBlockList[i].printElmtBlockInfo();
    }
    MessagePrinter::printStars();
}

void BulkElmtSystem::releaseMemory(){
    m_ElmtBlockList.clear();
    m_ElmtBlockNum=0;
    
    m_LocalElemental_ElmtBlockID.clear();
    m_GlobalBulkElmtsNum=0;
    m_LocalBulkElmtsNum=0;
}