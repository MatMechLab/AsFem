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
//+++ Date   : 2020.07.01
//+++ Purpose: Define the input block for the element in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iomanip>
#include <string>
#include <vector>

#include "ElmtSystem/ElmtType.h"
#include "MateSystem/MateType.h"

#include "Utils/MessagePrinter.h"

using namespace std;
class ElmtBlock{
public:
    ElmtBlock(){
        _DofsIDList.clear();
        _DofsNameList.clear();
        _nDofs=0;
        _ElmtBlockName.clear();
        _ElmtTypeName.clear();
        _MateBlockName.clear();
        _DomainName="alldomain";
        _ElmtType=ElmtType::NULLELMT;
        _MateType=MateType::NULLMATE;
        _MateIndex=0;
    }

    vector<int>    _DofsIDList;
    vector<string> _DofsNameList;
    int            _nDofs;
    string         _ElmtBlockName;
    string         _ElmtTypeName;
    string         _MateBlockName;
    string         _DomainName;
    ElmtType       _ElmtType=ElmtType::NULLELMT;
    MateType       _MateType=MateType::NULLMATE;
    int            _MateIndex=0;
    
    void Init(){
        _DofsIDList.clear();
        _DofsNameList.clear();
        _nDofs=0;
        _ElmtBlockName.clear();
        _ElmtTypeName.clear();
        _MateBlockName.clear();
        _DomainName="alldomain";
        _ElmtType=ElmtType::NULLELMT;
        _MateType=MateType::NULLMATE;
        _MateIndex=0;
    }

    void PrintInfo()const{
        // char buff[70];
        string str;
        MessagePrinter::PrintNormalTxt(" +sub bulk element block information summary:");
        str="   block name = "+_ElmtBlockName+", using ["+_MateBlockName+"] mate block ( index= "+to_string(_MateIndex)+" )";
        MessagePrinter::PrintNormalTxt(str);

        str="   DoFs name ( dofs= "+to_string(_nDofs)+" ) = ";
        for(auto it:_DofsNameList) str+=it+"  ";
        MessagePrinter::PrintNormalTxt(str);

        str="   DoFs ID ( dofs = "+to_string(_nDofs)+" ) = ";
        for(auto it:_DofsIDList) str+=to_string(it)+"  ";
        MessagePrinter::PrintNormalTxt(str);

        str="   domain name ="+_DomainName;
        MessagePrinter::PrintNormalTxt(str);
    }
};