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
//+++ Purpose: Implement the materials system for AsFem
//+++          this class offer some built-in material models
//+++          as well as the User-Defined-Material (umat) models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/MateSystem.h"

MateSystem::MateSystem(){
    _nMateBlocks=0;
    _MateBlockList.clear();
}

void MateSystem::AddMateBlock2List(MateBlock &mateblock){
    string msg;
    if(_MateBlockList.size()<1){
        _MateBlockList.push_back(mateblock);
        _nMateBlocks=1;
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_MateBlockList.size();i++){
            if(_MateBlockList[i]._MateBlockName==mateblock._MateBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            _MateBlockList.push_back(mateblock);
            _nMateBlocks=int(_MateBlockList.size());
        }
        else{
            msg="duplicated ["+mateblock._MateBlockName+"] in the [mates] sub block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}