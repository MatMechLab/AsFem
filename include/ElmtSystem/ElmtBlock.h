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

#include "ElmtSystem/ElmtType.h"

using namespace std;

using namespace std;
class ElmtBlock{
public:
    ElmtBlock(){
        _DofsIDList.clear();
        _DofsNameList.clear();
        _nDofs=0;
        _MateBlockName.clear();
        _DomainName="alldomain";
        _ElmtType=ElmtType::NULLELMT;
    }

    vector<int>    _DofsIDList;
    vector<string> _DofsNameList;
    int            _nDofs;
    string         _MateBlockName;
    string         _DomainName;
    ElmtType       _ElmtType=ElmtType::NULLELMT;
    
    void Init(){
        _DofsIDList.clear();
        _DofsNameList.clear();
        _nDofs.clear();
        _MateBlockName.clear();
        _DomainName="alldomain";
        _ElmtType=ElmtType::NULLELMT;
    }
};