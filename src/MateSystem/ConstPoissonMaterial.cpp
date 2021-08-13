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

void ConstPoissonMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}

}

//********************************************************************
void ConstPoissonMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}



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
    Mate.ScalarMaterials("sigma")=InputParams[0];// sigma
    Mate.ScalarMaterials("dsigmadu")=0.0;// dsigma/dphi
    Mate.ScalarMaterials("f")=InputParams[1];// F
    Mate.ScalarMaterials("dfdu")=0.0;// dF/dphi
    Mate.VectorMaterials("gradu")=elmtsoln.gpGradU[1];// the gradient of u

}
