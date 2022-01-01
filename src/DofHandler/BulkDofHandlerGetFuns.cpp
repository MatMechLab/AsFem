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
//+++ Date   : 2020.07.10
//+++ Purpose: Implement the getting functions for bulk dof handler
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

//********************************************************
int BulkDofHandler::GetDofIDviaDofName(string dofname)const{
    for(const auto &it:_DofName2IDList){
        if(it.first==dofname){
            return it.second;
        }
    }
    return -1;// for invalid dof id!
}
//********************************************************
vector<int> BulkDofHandler::GetDofsIndexFromNameVec(vector<string> namelist)const{
    vector<int> ids;
    ids.clear();
    if(IsValidDofNameVec(namelist)){
        for(const auto &it:namelist){
            for(const auto &jit:_DofName2IDList){
                if(it==jit.first){
                    ids.push_back(jit.second);
                }
            }
        }
    }
    else{
        MessagePrinter::PrintErrorTxt("invalid dofs namelist vec, you should use the name defined in [dofs] block");
        MessagePrinter::AsFem_Exit();
    }
    return ids;
}
