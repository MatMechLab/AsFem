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

void DoubleWellFreeEnergyMaterial::InitMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                          const Vector3d &gpCoord, const vector<double> &InputParams,
                                                          const vector<double> &gpU, const vector<double> &gpUdot,
                                                          const vector<Vector3d> &gpGradU,
                                                          const vector<Vector3d> &gpGradUdot, Materials &Mate) {
    // Here we do not consider any initial internal stats
    if(t||dt||nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){

    }
}

void DoubleWellFreeEnergyMaterial::ComputeF(const vector<double> &InputParams, const vector<double> &U,
                                            const vector<double> &dUdt, vector<double> &F) {
    if(dUdt[0]){}
    c=U[1];
    F[0]=c*log(c)+(1-c)*log(1-c)+InputParams[1]*c*(1-c);
}
void DoubleWellFreeEnergyMaterial::ComputedFdU(const vector<double> &InputParams, const vector<double> &U,
                                               const vector<double> &dUdt, vector<double> &dF) {
    if(dUdt[0]){}
    c=U[1];
    dF[0]=log(c)-log(1-c)+InputParams[1]*(1-2*c);
}
void DoubleWellFreeEnergyMaterial::Computed2FdU2(const vector<double> &InputParams, const vector<double> &U,
                                                 const vector<double> &dUdt, vector<double> &d2F) {
    if(dUdt[0]){}
    c=U[1];
    d2F[0]=1.0/c+1.0/(1-c)-2*InputParams[1];
}

void DoubleWellFreeEnergyMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                             const Vector3d &gpCoord, const vector<double> &InputParams,
                                                             const vector<double> &gpU, const vector<double> &gpUOld,
                                                             const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                                                             const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                                                             const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld,
                                                             const Materials &MateOld, Materials &Mate) {
    if(t||dt||nDim||gpCoord(1)||gpU[0]||gpUOld[0]||
       gpUdot[0]||gpUdotOld[0]||gpGradU[0](1)||gpGradUOld[0](1)||
       gpGradUdot[0](1)||gpGradUdotOld[0](1)||MateOld.ScalarMaterials.size()){}// get rid of unused warning

    if(InputParams.size()<3){
        MessagePrinter::PrintErrorTxt("for double well free energy material, three parameters are required, you need to give: D, Chi, and Kappa");
        MessagePrinter::AsFem_Exit();
    }

    ComputeF(InputParams,gpU,gpUdot,_F);
    ComputedFdU(InputParams,gpU,gpUdot,_dFdc);
    Computed2FdU2(InputParams,gpU,gpUdot,_d2Fdc2);

    c=gpU[1];
    Mate.ScalarMaterials["M"]=InputParams[0]*c*(1-c);   // M
    Mate.ScalarMaterials["dMdc"]=InputParams[0]*(1-2*c);// dM/dc

    Mate.ScalarMaterials["F"]=_F[0];
    Mate.ScalarMaterials["dFdc"]=_dFdc[0];
    Mate.ScalarMaterials["d2Fdc2"]=_d2Fdc2[0];

    Mate.ScalarMaterials["Kappa"]=InputParams[2];
}