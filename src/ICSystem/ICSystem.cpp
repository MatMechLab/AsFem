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
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) constant value ic
//+++               2) random value ic
//+++               3) other type or user defined ic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

ICSystem::ICSystem(){
    _nICBlocks=0;
    _ICBlockList.clear();
}

//*********************************************
void ICSystem::AddICBlock2List(ICBlock &icblock){
    string msg;
    if(_ICBlockList.size()<1){
        _ICBlockList.push_back(icblock);
        _nICBlocks=int(_ICBlockList.size());
    }
    else{
        bool NotInList=true;
        for(const auto &it:_ICBlockList){
            if(it._ICBlockName==icblock._ICBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            _ICBlockList.push_back(icblock);
            _nICBlocks=int(_ICBlockList.size());
        }
        else{
            msg="duplicated ["+icblock._ICBlockName+"] in the [ics] sub block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}
