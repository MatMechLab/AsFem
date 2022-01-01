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
//+++ Date   : 2021.04.10
//+++ Purpose: Calculate the material properties required by Diffusion
//+++          element. In this code, we can define:
//+++           1) D
//+++           2) dD/dc(=0)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/ConstDiffusionMaterial.h"

void ConstDiffusionMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}

}
//****************************************************************************
void ConstDiffusionMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}


    if(InputParams.size()<1){
        MessagePrinter::PrintErrorTxt("for constant diffusion material, one parameter, namely the diffusivity, is required");
        MessagePrinter::AsFem_Exit();
    }

    Mate.ScalarMaterials("D")=InputParams[0];// D
    Mate.ScalarMaterials("dDdc")=0.0;        // dD/dc

    Mate.VectorMaterials("gradc")=elmtsoln.gpGradU[1];

}

