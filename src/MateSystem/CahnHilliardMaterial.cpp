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
//+++ Date   : 2020.12.30
//+++ Purpose: Calculate the material properties required by Diffusion
//+++          element. In this code, we can define:
//+++           1) M
//+++           2) dM/dc
//+++           3) F
//+++           4) dF/dc
//+++           5) d2F/dc2
//+++           6) Kappa
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::CahnHilliardMaterial(const int &nDim, const double &t, const double &dt,
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
        MessagePrinter::PrintErrorTxt("for cahn hilliard material, three parameters are required, you need to give: D, Chi, and Kappa");
        MessagePrinter::AsFem_Exit();
    }

    //double c;
    //c=gpU[1];
    //_ScalarMaterials["M"]=InputParams[0]*c*(1-c);// M
    //_ScalarMaterials["dMdc"]=InputParams[0]*(1-2*c);// dM/dc
    //_ScalarMaterials["F"]=c*log(c)+(1-c)*log(1-c)+InputParams[1]*c*(1-c);
    //_ScalarMaterials["dFdc"]=log(c)-log(1-c)+InputParams[1]*(1-2*c);
    //_ScalarMaterials["d2Fdc2"]=1.0/c+1.0/(1-c)-2*InputParams[1];
    //_ScalarMaterials["Kappa"]=InputParams[2];
}

