//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.10
//+++ Purpose: Define the element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

BulkElmtSystem::BulkElmtSystem(){
    _nBulkElmtBlocks=0;
    _BulkElmtBlockList.clear();
}

void BulkElmtSystem::InitBulkElmtSystem(){
    _nBulkElmtBlocks=0;
    _BulkElmtBlockList.clear();
}

//***********************************
void BulkElmtSystem::AddBulkElmtBlock2List(ElmtBlock &elmtBlock){
    string msg;
    if(_BulkElmtBlockList.size()<1){
        _BulkElmtBlockList.push_back(elmtBlock);
        _nBulkElmtBlocks=1;
    }
    else{
        bool IsBlockNameUnique,IsTypeAndDomainUnique;
        IsBlockNameUnique=true;
        for(const auto &it:_BulkElmtBlockList){
            if(it._ElmtBlockName==elmtBlock._ElmtBlockName){
                IsBlockNameUnique=false;
                break;
            }
        }
        if(IsBlockNameUnique){
            // now the block name is unique, then we need to check the type and domain
            IsTypeAndDomainUnique=true;
            for(const auto &it:_BulkElmtBlockList){
                if(it._ElmtType==elmtBlock._ElmtType&&
                   it._DomainName==elmtBlock._DomainName){
                    IsTypeAndDomainUnique=false;
                    break;
                }
            }
            if(IsTypeAndDomainUnique){
                _BulkElmtBlockList.push_back(elmtBlock);
                _nBulkElmtBlocks+=1;
            }
            else{
                msg="duplicated elmt type and domain settings in the [elmts] block, the unique 'type=' and 'domain=' should be given";
                MessagePrinter::PrintErrorTxt(msg);
                MessagePrinter::AsFem_Exit();
            }
        }
        else{
            msg="duplicated ["+elmtBlock._ElmtBlockName+"] is found in [elmts] sub block, the block name should be unique";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}