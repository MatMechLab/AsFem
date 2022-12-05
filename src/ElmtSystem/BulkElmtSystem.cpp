//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.12
//+++ Purpose: the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

BulkElmtSystem::BulkElmtSystem(){
    m_elmtblock_num=0;
    m_elmtblock_list.clear();
    m_bulk_elmts=0;
    m_elemental_elmtblock_id.clear();
}

void BulkElmtSystem::addElmtBlock2List(ElmtBlock t_elmtblock){
    if(m_elmtblock_list.size()<1){
        m_elmtblock_list.push_back(t_elmtblock);
        m_elmtblock_num=1;
    }
    else{
        bool IsExist=false;
        for(int i=0;i<static_cast<int>(m_elmtblock_list.size());i++){
            if(t_elmtblock.m_elmt_blockname==m_elmtblock_list[i].m_elmt_blockname){
                IsExist=true;
                break;
            }
        }
        if(!IsExist){
            m_elmtblock_list.push_back(t_elmtblock);
            m_elmtblock_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt(t_elmtblock.m_elmt_blockname+" is already exist, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}

void BulkElmtSystem::printBulkElmtSystemInfo()const{
    MessagePrinter::printNormalTxt("Bulk element system information summary");
    MessagePrinter::printNormalTxt(" number of bulk element blocks = "+to_string(getBulkElmtBlocksNum()));
    for(int i=0;i<m_elmtblock_num;i++){
        m_elmtblock_list[i].printElmtBlockInfo();
    }
    MessagePrinter::printStars();
}

void BulkElmtSystem::releaseMemory(){
    m_elmtblock_list.clear();
    m_elmtblock_num=0;
    
    m_elemental_elmtblock_id.clear();
    m_bulk_elmts=0;
}