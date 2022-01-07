//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
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
    
    ca=InputParams[1];
    cb=InputParams[2];
    factor=InputParams[3];
    
    F[0]=factor*(c-ca)*(c-ca)*(c-cb)*(c-cb);
}
//****************************************************************

void DoubleWellFreeEnergyMaterial::ComputedFdU(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &dF){
    
    c=elmtsoln.gpU[1];
    
    ca=InputParams[1];
    cb=InputParams[2];
    factor=InputParams[3];

    dF[0]=factor*2*(c-ca)*(c-cb)*(2*c-ca-cb);
}

void DoubleWellFreeEnergyMaterial::Computed2FdU2(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &d2F){
    c=elmtsoln.gpU[1];
    
    ca=InputParams[1];
    cb=InputParams[2];
    factor=InputParams[3];

    d2F[0]=2*(ca*ca+4*ca*cb-6*ca*c+cb*cb-6*cb*c+6*c*c);
}

void DoubleWellFreeEnergyMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //********************************************
    // get rid of unused warnings
    //********************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}


    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for double well free energy material(f=h(c-ca)^2(c-cb)^2), five parameters are required, you need to give: D, ca, cb, h,  and Kappa");
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

    Mate.ScalarMaterials("Kappa")=InputParams[4];

    Mate.VectorMaterials("gradc") =elmtsoln.gpGradU[1];// for output
    Mate.VectorMaterials("gradmu")=elmtsoln.gpGradU[2];// for output

}
