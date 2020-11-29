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

#include "ElmtSystem/ElmtBlock.h"
#include "ElmtSystem/ElmtType.h"
#include "FESystem/FECalcType.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;

class BulkElmtSystem{
public:
    BulkElmtSystem();

    void InitBulkElmtSystem();


    void AddBulkElmtBlock2List(ElmtBlock &elmtBlock);
    ElmtBlock GetIthBulkElmtBlock(const int &i)const{
        return _BulkElmtBlockList[i-1];
    }
    inline int GetBulkElmtBlockNums()const{
        return _nBulkElmtBlocks;
    }


    void RunBulkElmtLibs(const FECalcType &calctype,const ElmtType &elmtytype,
                        const int &nDim,const int &nNodes,
                        const double &t,const double &dt,const double (&ctan)[2],
                        const Vector3d &gpCoords,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                        // const ShapeFun &shp,
                        const vector<double> &ScalarMaterials,
                        const vector<Vector3d> &VectorMaterials,
                        const vector<RankTwoTensor> &Rank2Materials,
                        const vector<RankFourTensor> &Rank4Materials,
                        vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                        MatrixXd &localK,VectorXd &localR);

protected:
    int _nBulkElmtBlocks;
    vector<ElmtBlock> _BulkElmtBlockList;
};