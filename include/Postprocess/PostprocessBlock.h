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
    }

    int _NodeID;
    int _ElementID;
    vector<string> _DomainNameList;
    vector<string> _BoundaryNameList;
    string _PPSBlockName;
    PostprocessType _PostprocessType;
    string _PostprocessTypeName;
    string _VariableName;// VariableName could be name in [dofs] block, or name ini [projection] block, or even our material properties

    void Init(){
        _NodeID=-1;
        _ElementID=-1;
        _DomainNameList.clear();
        _BoundaryNameList.clear();

        _PPSBlockName.clear();

        _PostprocessType=PostprocessType::NULLPPS;
        _PostprocessTypeName="none";
        _VariableName.clear();
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
        if(_NodeID>0){
            str="   node id="+to_string(_NodeID);
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_ElementID>0){
            str="   element id="+to_string(_NodeID);
            MessagePrinter::PrintNormalTxt(str);
        }
    }
};