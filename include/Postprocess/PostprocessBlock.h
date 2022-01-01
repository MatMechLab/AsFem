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
//+++ Date   : 2021.02.21
//+++ Purpose: define the postprocess block for input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "Postprocess/PostprocessType.h"
#include "Utils/MessagePrinter.h"

using namespace std;

class PostprocessBlock{
public:
    PostprocessBlock(){
        _NodeID=-1;
        _ElementID=-1;
        _DomainNameList.clear();
        _BoundaryNameList.clear();

        _PPSBlockName.clear();

        _PostprocessType=PostprocessType::NULLPPS;
        _PostprocessTypeName="none";
        _VariableName.clear();

        _ProjVariableName.clear();
        _ScalarMateName.clear();_VectorMateName.clear();_Rank2MateName.clear();_Rank4MateName.clear();
        _iInd=-1;_jInd=-1;_Component=-1;
    }

    int _NodeID;
    int _ElementID;
    vector<string> _DomainNameList;
    vector<string> _BoundaryNameList;
    string _PPSBlockName;
    PostprocessType _PostprocessType;
    string _PostprocessTypeName;
    string _VariableName;// VariableName could be name in [dofs] block, or name ini [projection] block, or even our material properties
    string _ProjVariableName;
    string _ScalarMateName,_VectorMateName,_Rank2MateName,_Rank4MateName;
    int _iInd,_jInd,_Component;

    void Init(){
        _NodeID=-1;
        _ElementID=-1;
        _DomainNameList.clear();
        _BoundaryNameList.clear();

        _PPSBlockName.clear();

        _PostprocessType=PostprocessType::NULLPPS;
        _PostprocessTypeName="none";
        _VariableName.clear();

        _ProjVariableName.clear();
        _ScalarMateName.clear();_VectorMateName.clear();_Rank2MateName.clear();_Rank4MateName.clear();
        _iInd=-1;_jInd=-1;_Component=-1;
    }

    void PrintInfo()const{
        string str;
        MessagePrinter::PrintNormalTxt(" +sub postprocess block information summary:");

        str="   block name = "+_PPSBlockName+", pps type="+_PostprocessTypeName;
        MessagePrinter::PrintNormalTxt(str);
        if(_DomainNameList.size()>0){
            str="   domain=";
            for(auto it:_DomainNameList) str+=it+" ";
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_BoundaryNameList.size()>0){
            str="   boundary=";
            for(auto it:_BoundaryNameList) str+=it+" ";
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_VariableName.size()>0){
            str="   variable="+_VariableName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_ProjVariableName.size()>0){
            str="   projected variable="+_ProjVariableName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_ScalarMateName.size()>0){
            str="   scalar material="+_ScalarMateName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_VectorMateName.size()>0){
            str="   vector material="+_VectorMateName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_Rank2MateName.size()>0){
            str="   rank-2 tensor material="+_Rank2MateName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_Rank4MateName.size()>0){
            str="   rank-4 tensor material="+_Rank4MateName;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_NodeID>0){
            str="   node id="+to_string(_NodeID);
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_ElementID>0){
            str="   element id="+to_string(_NodeID);
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_iInd>0){
            str="   i-index ="+to_string(_iInd);
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_jInd>0){
            str="   j-index ="+to_string(_jInd);
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_Component>0){
            str="   component ="+to_string(_Component);
            MessagePrinter::PrintNormalTxt(str);
        }
    }
};