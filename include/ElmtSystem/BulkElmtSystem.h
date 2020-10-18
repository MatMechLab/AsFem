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
//+++ Date   : 2020.10.18
//+++ Purpose: Define the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include <iostream>
#include <iomanip>

#include "Utils/MessagePrinter.h"

#include "ElmtSystem/ElmtType.h"
#include "ElmtSystem/ElmtBlock.h"

using namespace std;

class BulkElmtSystem{
public:
    BulkElmtSystem();


    enum BulkElmtCalcType{
        ComputeResidual,
        ComputeJacobian,
        InitHistoryValue,
        UpdateHistoryValue,
        ComputeProjectionValue
    };

    void AddElmtBlock2List(ElmtBlock &elmtBlock);
    ElmtBlock GetIthElmtBlock(const int &i)const{
        return _ElmtBlockList[i-1];
    }

protected:
    int _nElmtBlocks;
    vector<ElmtBlock> _ElmtBlockList;
};