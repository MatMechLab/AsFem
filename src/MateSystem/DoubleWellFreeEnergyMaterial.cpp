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
//+++ Date   : 2021.04.04
//+++ Purpose: Calculate the free energy, chemical potential and its
//+++          derivatives of double well free energy material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/DoubleWellFreeEnergyMaterial.h"

DoubleWellFreeEnergyMaterial::DoubleWellFreeEnergyMaterial() {
    _F.resize(1,0.0);
    _dFdc.resize(1,0.0);
    _d2Fdc2.resize(1,0.0);
}
//****************************************************************
void DoubleWellFreeEnergyMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal stats
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()||Mate.GetScalarMate().size()){}


}
//*****************************************************************
void DoubleWellFreeEnergyMaterial::ComputeF(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &F){
    c=elmtsoln.gpU[1];
    F[0]=c*log(c)+(1-c)*log(1-c)+InputParams[1]*c*(1-c);
}
//****************************************************************

void DoubleWellFreeEnergyMaterial::ComputedFdU(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &dF){
    
    c=elmtsoln.gpU[1];
    dF[0]=log(c)-log(1-c)+InputParams[1]*(1-2*c);
}

void DoubleWellFreeEnergyMaterial::Computed2FdU2(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &d2F){
    c=elmtsoln.gpU[1];
    d2F[0]=1.0/c+1.0/(1-c)-2*InputParams[1];
}

void DoubleWellFreeEnergyMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //********************************************
    // get rid of unused warnings
    //********************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}


    if(InputParams.size()<3){
        MessagePrinter::PrintErrorTxt("for double well free energy material, three parameters are required, you need to give: D, Chi, and Kappa");
        MessagePrinter::AsFem_Exit();
    }

    ComputeF(InputParams,elmtsoln,_F);
    ComputedFdU(InputParams,elmtsoln,_dFdc);
    Computed2FdU2(InputParams,elmtsoln,_d2Fdc2);

    c=elmtsoln.gpU[1];
    Mate.ScalarMaterials("M")=InputParams[0]*c*(1-c);   // M
    Mate.ScalarMaterials("dMdc")=InputParams[0]*(1-2*c);// dM/dc

    Mate.ScalarMaterials("F")=_F[0];
    Mate.ScalarMaterials("dFdc")=_dFdc[0];
    Mate.ScalarMaterials("d2Fdc2")=_d2Fdc2[0];

    Mate.ScalarMaterials("Kappa")=InputParams[2];

    Mate.VectorMaterials("gradc") =elmtsoln.gpGradU[1];// for output
    Mate.VectorMaterials("gradmu")=elmtsoln.gpGradU[2];// for output

}
