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
//+++ Purpose: Calculate the material properties required by Diffusion
//+++          element. In this code, we can define:
//+++           1) D
//+++           2) dD/dc(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/ConstDiffusionMaterial.h"

void ConstDiffusionMaterial::InitMaterialProperties(const int &nDim, const Vector3d &gpCoord,
                                                    const vector<double> &InputParams, const vector<double> &gpU,
                                                    const vector<double> &gpUdot, const vector<Vector3d> &gpGradU,
                                                    const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){}
}
//****************************************************************************
void ConstDiffusionMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
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

    if(InputParams.size()<1){
        MessagePrinter::PrintErrorTxt("for constant diffusion material, one parameter, namely the diffusivity, is required");
        MessagePrinter::AsFem_Exit();
    }

    Mate.ScalarMaterials["D"]=InputParams[0];// D
    Mate.ScalarMaterials["dDdc"]=0.0;        // dD/dc

    Mate.VectorMaterials["gradc"]=gpGradU[1];
}

