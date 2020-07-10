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

#include "ElmtSystem/ElmtSystem.h"

ElmtSystem::ElmtSystem(){
    _nElmtBlocks=0;
    _ElmtBlockList.clear();
}

//***********************************
void ElmtSystem::AddElmtBlock2List(ElmtBlock &elmtBlock){
    string msg;
    if(_ElmtBlockList.size()<1){
        _ElmtBlockList.push_back(elmtBlock);
    }
    else{
        bool IsBlockNameUnique,IsTypeAndDomainUnique;
        IsBlockNameUnique=true;
        for(const auto &it:_ElmtBlockList){
            if(it._ElmtBlockName==elmtBlock._ElmtBlockName){
                IsBlockNameUnique=false;
                break;
            }
        }
        if(IsBlockNameUnique){
            // now the block name is unique, then we need to check the type and domain
            IsTypeAndDomainUnique=true;
            for(const auto &it:_ElmtBlockList){
                if(it._ElmtType==elmtBlock._ElmtType&&
                   it._DomainName==elmtBlock._DomainName){
                    IsTypeAndDomainUnique=false;
                    break;
                }
            }
            if(IsTypeAndDomainUnique){
                _ElmtBlockList.push_back(elmtBlock);
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