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
//+++ Date   : 2020.07.12
//+++ Purpose: define the solution system for AsFem, where all the 
//+++          results/material variables/projection quantities 
//+++          should be stored here by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

SolutionSystem::SolutionSystem(){

    _IsInit=false;_IsProjection=false;
    _DofNameList.clear();
    _ProjectionNameList.clear();

    _nHistPerGPoint=6;_nProjPerNode=9;
    _nGPointsPerBulkElmt=0;
    _nDofs=0;_nNodes=0;_nElmts=0;
    _HasProjNameList=false;
    _HasDofNameList=false;
}

void SolutionSystem::PrintProjectionInfo()const{
    if(_ProjectionNameList.size()<1) return;
    MessagePrinter::PrintNormalTxt("Projection information summary:");
    string msg;
    msg="projected variables: ";
    for(const auto it:_ProjectionNameList){
        msg+=it+" ";
    }
    MessagePrinter::PrintNormalTxt(msg);
    MessagePrinter::PrintDashLine();
}