//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
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

    _nHistPerGPoint=6;
    _nProjPerNode=0;
    _nScalarProjPerNode=0;_nVectorProjPerNode=0;_nRank2ProjPerNode=0;_nRank4ProjPerNode=0;
    _nGPointsPerBulkElmt=0;
    _nDofs=0;_nNodes=0;_nElmts=0;
    _HasProjNameList=false;
    _HasDofNameList=false;

    //***************************************
    _IsInit=false;_IsProjection=false;
    _DofNameList.clear();
    _ProjectionNameList.clear();
    _ScalarMateProjectionNameList.clear();_VectorMateProjctionNameList.clear();
    _Rank2MateProjectionNameList.clear();_Rank4MateProjectionNameList.clear();

    _nHistPerGPoint=6;_nProjPerNode=0;
    _nScalarProjPerNode=0;_nVectorProjPerNode=0;_nRank2ProjPerNode=0;_nRank4ProjPerNode=0;
    _nGPointsPerBulkElmt=0;
    _nDofs=0;_nNodes=0;_nElmts=0;
    _HasProjNameList=false;
    _HasScalarMateProjName=false;
    _HasVectorMateProjName=false;
    _HasRank2MateProjName=false;
    _HasRank4MateProjName=false;
    _HasDofNameList=false;
}

void SolutionSystem::PrintProjectionInfo()const{
    MessagePrinter::PrintNormalTxt("Projection information summary:");
    string msg;
    msg="projected variables: ";
    for(const auto &it:_ProjectionNameList){
        msg+=it+" ";
    }
    if(GetProjNumPerNode()>0) MessagePrinter::PrintNormalTxt(msg);
    //********************************************
    msg="projected scalar material: ";
    for(const auto &it:_ScalarMateProjectionNameList){
        msg+=it+" ";
    }
    if(GetScalarMateProjNumPerNode()>0) MessagePrinter::PrintNormalTxt(msg);

    //********************************************
    msg="projected vector material: ";
    for(const auto &it:_VectorMateProjctionNameList){
        msg+=it+" ";
    }
    if(GetVectorMateProjNumPerNode()>0) MessagePrinter::PrintNormalTxt(msg);

    //********************************************
    msg="projected rank-2 tensor material: ";
    for(const auto &it:_Rank2MateProjectionNameList){
        msg+=it+" ";
    }
    if(GetRank2MateProjNumPerNode()>0) MessagePrinter::PrintNormalTxt(msg);

    //********************************************
    msg="projected rank-4 tensor material: ";
    for(const auto &it:_Rank4MateProjectionNameList){
        msg+=it+" ";
    }
    if(GetRank4MateProjNumPerNode()>0) MessagePrinter::PrintNormalTxt(msg);

    MessagePrinter::PrintDashLine();
}

//*******************************************
void SolutionSystem::ReleaseMem(){
    VecDestroy(&_Unew);
    VecDestroy(&_U);
    VecDestroy(&_Utemp);
    VecDestroy(&_V);

    VecDestroy(&_Uold);
    VecDestroy(&_Vold);

    VecDestroy(&_Proj);

    VecDestroy(&_ProjScalarMate);
    VecDestroy(&_ProjVectorMate);
    VecDestroy(&_ProjRank2Mate);
    VecDestroy(&_ProjRank4Mate);

}
