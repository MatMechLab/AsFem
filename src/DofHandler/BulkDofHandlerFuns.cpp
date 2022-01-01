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
//+++ Purpose: Implement some basic check functions for dofhandler
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

//***************************************************
bool BulkDofHandler::IsValidDofName(string dofname)const{
    bool IsValid=false;
    for(const auto &it:_DofNameList){
        if(it==dofname){
            IsValid=true;
            break;
        }
    }
    return IsValid;
}
//***************************************************
bool BulkDofHandler::IsValidDofNameVec(vector<string> namelist)const{
    bool IsValid=true;
    for(const auto &it:namelist){
        IsValid=false;
        for(const auto &itj:_DofNameList){
            if(it==itj){
                IsValid=true;
                break;
            }
        }
        if(!IsValid){
            return false;
        }
    }
    return true;
}