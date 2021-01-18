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
//+++ Date   : 2021.01.18
//+++ Purpose: Calculate the User-Defined-Material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::User1Material(const int &nDim, const double &t, const double &dt,
                                   const vector<double> &InputParams, const Vector3d &gpCoord,
                                   const vector<double> &gpU, const vector<double> &gpV,
                                   const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradV,
                                   vector<double> &gpHist, const vector<double> &gpHistOld) {

    //*****************************************************************************
    //*** just to get rid of warnings, normal users dont need to do this
    //*****************************************************************************
    if(nDim||t||dt||InputParams.size()||
       gpCoord(1)||gpU.size()||gpV.size()||gpGradU.size()||gpGradV.size()||
       gpHist.size()||gpHistOld.size()){}

    if(InputParams.size()<3){
        MessagePrinter::PrintErrorTxt("for linear elastic material, two parameters are required, you need to give: E, nu, and D");
        MessagePrinter::AsFem_Exit();
    }

    const double E=InputParams[0];
    const double nu=InputParams[1];
    const double D=InputParams[2];

    // for our elasticity tensor
    _Rank4Materials["elasticity_tensor"].SetToZeros();
    _Rank4Materials["elasticity_tensor"].SetFromEandNu(E,nu);

    // for our small strain
    _Rank2Materials["strain"].SetToZeros();
    RankTwoTensor GradU;
    GradU.SetToZeros();
    if(nDim==1){
        GradU.SetFromGradU(gpGradU[2]);
    }
    else if(nDim==2){
        GradU.SetFromGradU(gpGradU[2],gpGradU[3]);
    }
    else if(nDim==3){
        GradU.SetFromGradU(gpGradU[2],gpGradU[3],gpGradU[4]);
    }
    // epsilon_ij=0.5*(Ui,j+Uj,i)
    _Rank2Materials["strain"]=(GradU+GradU.Transpose())*0.5;

    // our stress
    _Rank2Materials["stress"]=_Rank4Materials["elasticity_tensor"].DoubleDot(_Rank2Materials["strain"]);

    // for vonMises stress
    RankTwoTensor I,devStress;
    double trace;
    I.SetToIdentity();
    trace=_Rank2Materials["stress"].Trace();

    devStress=_Rank2Materials["stress"]-I*(trace/3.0);
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    _ScalarMaterials["vonMises"]=sqrt(1.5*devStress.DoubleDot(devStress));
    _ScalarMaterials["D"]=D;
}

