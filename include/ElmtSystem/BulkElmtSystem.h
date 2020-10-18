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

    void InitBulkElmtSystem();

    enum BulkElmtCalcType{
        ComputeResidual,
        ComputeJacobian,
        InitHistoryValue,
        UpdateHistoryValue,
        ComputeProjectionValue
    };

    void AddBulkElmtBlock2List(ElmtBlock &elmtBlock);
    ElmtBlock GetIthBulkElmtBlock(const int &i)const{
        return _BulkElmtBlockList[i-1];
    }
    inline int GetBulkElmtBlockNums()const{
        return _nBulkElmtBlocks;
    }

protected:
    int _nBulkElmtBlocks;
    vector<ElmtBlock> _BulkElmtBlockList;
};