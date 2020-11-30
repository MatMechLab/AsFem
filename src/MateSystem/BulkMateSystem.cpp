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
//+++ Date   : 2020.11.30
//+++ Purpose: Implement the materials system for our bulk element
//+++          in AsFem. It is different from another one, namely
//+++          the interface material system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

BulkMateSystem::BulkMateSystem(){
    _nBulkMateBlocks=0;
    _BulkMateBlockList.clear();
    _ScalarMaterials.clear();
    _VectorMaterials.clear();
    _Rank2Materials.clear();
    _Rank4Materials.clear();
}

//***************************************************
void BulkMateSystem::AddBulkMateBlock2List(MateBlock &mateblock){
    string msg;
    if(_BulkMateBlockList.size()<1){
        _BulkMateBlockList.push_back(mateblock);
        _nBulkMateBlocks=1;
    }
    else{
        bool NotInList=true;
        for(unsigned int i=0;i<_BulkMateBlockList.size();i++){
            if(_BulkMateBlockList[i]._MateBlockName==mateblock._MateBlockName){
                NotInList=false;
                break;
            }
        }
        if(NotInList){
            _BulkMateBlockList.push_back(mateblock);
            _nBulkMateBlocks=static_cast<int>(_BulkMateBlockList.size());
        }
        else{
            msg="duplicated ["+mateblock._MateBlockName+"] in the [mates] sub block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}
//*********************************************************
void BulkMateSystem::InitBulkMateSystem(){
    _ScalarMaterials.clear();
    _VectorMaterials.clear();
    _Rank2Materials.clear();
    _Rank4Materials.clear();
}
//***********************************************************
void BulkMateSystem::PrintBulkMateSystemInfo()const{
    if(_nBulkMateBlocks>0){
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        MessagePrinter::PrintDashLine();
        MessagePrinter::PrintNormalTxt("Material system information summary:");
        string str;
        for(auto it:_BulkMateBlockList){
            str="  mate block name="+it._MateBlockName+", mate type="+it._MateTypeName;
            MessagePrinter::PrintNormalTxt(str);
            str="  parameters=";
            for(auto sit:it._Parameters) str+=to_string(sit);
            MessagePrinter::PrintNormalTxt(str);
        }
    }
}