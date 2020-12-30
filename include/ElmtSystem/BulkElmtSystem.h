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
#include "MateSystem/MateSystem.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;

class MateSystem;

class BulkElmtSystem{
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
                        const int &nDim,const int &nNodes,
                        const double &t,const double &dt,const double (&ctan)[2],
                        const Vector3d &gpCoords,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                        const double &test,const double &trial,
                        const Vector3d &grad_test,const Vector3d &grad_trial,
                        const ScalarMateType &ScalarMaterials,
                        const VectorMateType &VectorMaterials,
                        const Rank2MateType &Rank2Materials,
                        const Rank4MateType &Rank4Materials,
                        vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                        MatrixXd &localK,VectorXd &localR);

    void PrintBulkElmtInfo()const;

protected:
    int _nBulkElmtBlocks;
    vector<ElmtBlock> _BulkElmtBlockList;


protected:
    //****************************************************************************
    //*** For AsFem's built-in elements and User-Defined-Elements (UEL)
    //*** One can implement his own model here
    //****************************************************************************
    //*** the laplace means: \int(\nabla U*\nabla test)dV
    void LaplaceElmt(const FECalcType &calctype,
                const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoords,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const double &test,const double &trial,
                const Vector3d &grad_test,const Vector3d &grad_trial,
                const ScalarMateType &ScalarMaterials,
                const VectorMateType &VectorMaterials,
                const Rank2MateType &Rank2Materials,
                const Rank4MateType &Rank4Materials,
                vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //*** the body source means  int(f*test)dV
    void BodySourceElmt(const FECalcType &calctype,
                const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoords,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const double &test,const double &trial,
                const Vector3d &grad_test,const Vector3d &grad_trial,
                const ScalarMateType &ScalarMaterials,
                const VectorMateType &VectorMaterials,
                const Rank2MateType &Rank2Materials,
                const Rank4MateType &Rank4Materials,
                vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //************************************************************************************
    //*** for general time derivative
    //************************************************************************************
    void TimeDerivElmt(const FECalcType &calctype,
                     const int &nDim,const int &nNodes,
                     const double &t,const double &dt,const double (&ctan)[2],
                     const Vector3d &gpCoords,
                     const vector<double> &gpU,const vector<double> &gpV,
                     const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                     const double &test,const double &trial,
                     const Vector3d &grad_test,const Vector3d &grad_trial,
                     const ScalarMateType &ScalarMaterials,
                     const VectorMateType &VectorMaterials,
                     const Rank2MateType &Rank2Materials,
                     const Rank4MateType &Rank4Materials,
                     vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                     MatrixXd &localK,VectorXd &localR);
    //************************************************************************************
    //*** for general poisson element
    //************************************************************************************
    void PoissonElmt(const FECalcType &calctype,
                const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoords,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const double &test,const double &trial,
                const Vector3d &grad_test,const Vector3d &grad_trial,
                const ScalarMateType &ScalarMaterials,
                const VectorMateType &VectorMaterials,
                const Rank2MateType &Rank2Materials,
                const Rank4MateType &Rank4Materials,
                vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);

    //************************************************************************************
    //*** for general diffusion element
    //************************************************************************************
    void DiffusionElmt(const FECalcType &calctype,
                     const int &nDim,const int &nNodes,
                     const double &t,const double &dt,const double (&ctan)[2],
                     const Vector3d &gpCoords,
                     const vector<double> &gpU,const vector<double> &gpV,
                     const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                     const double &test,const double &trial,
                     const Vector3d &grad_test,const Vector3d &grad_trial,
                     const ScalarMateType &ScalarMaterials,
                     const VectorMateType &VectorMaterials,
                     const Rank2MateType &Rank2Materials,
                     const Rank4MateType &Rank4Materials,
                     vector<double> &gpHist,vector<double> &gpHistOld,vector<double> &gpProj,
                     MatrixXd &localK,VectorXd &localR);


};