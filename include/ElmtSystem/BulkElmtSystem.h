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
#include "ElmtSystem/KobayashiElmt.h"
#include "ElmtSystem/StressDiffusionElmt.h"
#include "ElmtSystem/DiffusionFractureElmt.h"
#include "ElmtSystem/MechanicsCahnHilliardElmt.h"
#include "ElmtSystem/AllenCahnFractureElmt.h"
// for UEL
#include "ElmtSystem/User1Elmt.h"
#include "ElmtSystem/User2Elmt.h"
#include "ElmtSystem/User3Elmt.h"
#include "ElmtSystem/User4Elmt.h"
#include "ElmtSystem/User5Elmt.h"
#include "ElmtSystem/User6Elmt.h"
#include "ElmtSystem/User7Elmt.h"
#include "ElmtSystem/User8Elmt.h"
#include "ElmtSystem/User9Elmt.h"
#include "ElmtSystem/User10Elmt.h"
#include "ElmtSystem/User11Elmt.h"
#include "ElmtSystem/User12Elmt.h"
#include "ElmtSystem/User13Elmt.h"
#include "ElmtSystem/User14Elmt.h"
#include "ElmtSystem/User15Elmt.h"
#include "ElmtSystem/User16Elmt.h"
#include "ElmtSystem/User17Elmt.h"
#include "ElmtSystem/User18Elmt.h"
#include "ElmtSystem/User19Elmt.h"
#include "ElmtSystem/User20Elmt.h"

using namespace std;

class MateSystem;

class BulkElmtSystem: public PoissonElmt,
                      public DiffusionElmt,
                      public MechanicsElmt,
                      public CahnHilliardElmt,
                      public MechanicsCahnHilliardElmt,
                      public MieheFractureElmt,
                      public KobayashiElmt,
                      public StressDiffusionElmt,
                      public DiffusionFractureElmt,
                      public AllenCahnFractureElmt,
                      public User1Elmt,
                      public User2Elmt,
                      public User3Elmt,
                      public User4Elmt,
                      public User5Elmt,
                      public User6Elmt,
                      public User7Elmt,
                      public User8Elmt,
                      public User9Elmt,
                      public User10Elmt,
                      public User11Elmt,
                      public User12Elmt,
                      public User13Elmt,
                      public User14Elmt,
                      public User15Elmt,
                      public User16Elmt,
                      public User17Elmt,
                      public User18Elmt,
                      public User19Elmt,
                      public User20Elmt{
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
                         const double (&ctan)[3],
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
