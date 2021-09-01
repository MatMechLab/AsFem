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
//+++ Date   : 2020.10.18
//+++ Purpose: Define the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include <iostream>
#include <iomanip>

#include "Utils/MessagePrinter.h"

#include "ElmtSystem/ElmtBlock.h"
#include "ElmtSystem/ElmtType.h"
#include "FESystem/FECalcType.h"
#include "MateSystem/MateSystem.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

// for all the user-defined-elements
#include "ElmtSystem/PoissonElmt.h"
#include "ElmtSystem/DiffusionElmt.h"
#include "ElmtSystem/MechanicsElmt.h"
#include "ElmtSystem/CahnHilliardElmt.h"
#include "ElmtSystem/MieheFractureElmt.h"

using namespace std;

class MateSystem;

class BulkElmtSystem: public PoissonElmt,
                      public DiffusionElmt,
                      public MechanicsElmt,
                      public CahnHilliardElmt,
                      public MieheFractureElmt{
public:
    BulkElmtSystem();

    void InitBulkElmtSystem();
    void InitBulkElmtMateInfo(MateSystem &matesystem);

    void AddBulkElmtBlock2List(ElmtBlock &elmtBlock);
    ElmtBlock GetIthBulkElmtBlock(const int &i)const{
        return _BulkElmtBlockList[i-1];
    }
    inline int GetBulkElmtBlockNums()const{
        return _nBulkElmtBlocks;
    }


    void RunBulkElmtLibs(const FECalcType &calctype,const ElmtType &elmtytype,
                         const double (&ctan)[2],
                         const LocalElmtInfo &elmtinfo,
                         const LocalElmtSolution &soln,
                         const LocalShapeFun &shp,
                         const Materials &Mate,const Materials &MateOld,
                         ScalarMateType &gpProj,
                         MatrixXd &localK,VectorXd &localR);

    void PrintBulkElmtInfo()const;

protected:
    int _nBulkElmtBlocks;
    vector<ElmtBlock> _BulkElmtBlockList;

protected:

};
