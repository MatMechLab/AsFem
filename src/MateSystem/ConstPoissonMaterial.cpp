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
//+++ Date   : 2021.04.10
//+++ Purpose: Calculate the material properties required by Poisson
//+++          element. In this code, we can define:
//+++           1) Sigma
//+++           2) dSigma/du(=0)
//+++           3) F
//+++           4) dF/du(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/ConstPoissonMaterial.h"

void ConstPoissonMaterial::InitMaterialProperties(const int &nDim, const Vector3d &gpCoord,
                                                  const vector<double> &InputParams, const vector<double> &gpU,
                                                  const vector<double> &gpUdot, const vector<Vector3d> &gpGradU,
                                                  const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
    gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){}
}

void ConstPoissonMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                     const Vector3d &gpCoord, const vector<double> &InputParams,
                                                     const vector<double> &gpU, const vector<double> &gpUOld,
                                                     const vector<double> &gpUdot, const vector<double> &gpUdotOld,
                                                     const vector<Vector3d> &gpGradU,
                                                     const vector<Vector3d> &gpGradUOld,
                                                     const vector<Vector3d> &gpGradUdot,
                                                     const vector<Vector3d> &gpGradUdotOld, const Materials &MateOld,
                                                     Materials &Mate) {
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(t||dt||nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUOld[0]||gpUdot[0]||gpUdotOld[0]||
       gpGradU[0](1)||gpGradUOld[0](1)||gpGradUdot[0](1)||gpGradUdotOld[0](1)||
       MateOld.ScalarMaterials.size()||Mate.ScalarMaterials.size()){}

    if(InputParams.size()<2){
        MessagePrinter::PrintErrorTxt("for const poisson material, two parameters are required. sigma*div(grad(phi))=F, so sigma and F are required");
        MessagePrinter::AsFem_Exit();
    }

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    //**** MateVals[0]-->store sigma
    //**** MateVals[1]-->store dsigma/dphi(for constant case, it is zero)
    //**** MateVals[2]-->store F
    //**** MateVals[3]-->store dF/dphi (for constant case, it is zero)
    Mate.ScalarMaterials["sigma"]=InputParams[0];// sigma
    Mate.ScalarMaterials["dsigmadu"]=0.0;// dsigma/dphi
    Mate.ScalarMaterials["f"]=InputParams[1];// F
    Mate.ScalarMaterials["dfdu"]=0.0;// dF/dphi
    Mate.VectorMaterials["gradu"]=gpGradU[1];
}