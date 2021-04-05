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
//+++           1) D
//+++           2) dD/dc(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::ConstDiffusionMaterial(const int &nDim, const double &t, const double &dt,
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

    if(InputParams.size()<1){
        MessagePrinter::PrintErrorTxt("for const diffusion material, one parameter is required, the diffusivity is required");
        MessagePrinter::AsFem_Exit();
    }

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    //**** MateVals[0]-->store D
    //**** MateVals[1]-->store dD/dc(for constant case, it is zero)
    //_ScalarMaterials["D"]=InputParams[0];// sigma
    //_ScalarMaterials["dDdc"]=0.0;// dsigma/dphi
}

