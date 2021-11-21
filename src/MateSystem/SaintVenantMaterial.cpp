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
//+++ Date   : 2021.10.31
//+++ Purpose: Implement the calculation of SaintVenant 
//+++          hyperelastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/SaintVenantMaterial.h"


void SaintVenantMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();
    Mate.Rank2Materials("F").SetToIdentity();
    Mate.Rank2Materials("PK1").SetToZeros();
    Mate.Rank2Materials("PK2").SetToZeros();
}

//*********************************************************
void SaintVenantMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain){
    // here we calculate the deformation gradient tensor as well as euler-lagrange strain tensor
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1]);
    }
    else if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    _I.SetToIdentity();
    _F=_GradU+_I;// F=I+U_{i,j}
    _C=_F.Transpose()*_F; //C=F^T*F
    Strain=(_C-_I)*0.5;

}
//****************************************************************
void SaintVenantMaterial::ComputeStressAndJacobian(const vector<double> &InputParams, const RankTwoTensor &Strain, RankTwoTensor &Stress, RankFourTensor &Jacobian){
    double E=InputParams[0];
    double nu=InputParams[1];

    Jacobian.SetFromEandNu(E,nu);
    _pk2=Jacobian.DoubleDot(Strain);
    Stress=_F*_pk2;// convert it to 1st Piola-Kirchhoff stress

}
//***************************************************************
void SaintVenantMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    if(InputParams.size()<2){
        MessagePrinter::PrintErrorTxt("for the SaintVenantMaterial, two parameters are required, you need to give: E and nu");
        MessagePrinter::AsFem_Exit();
    }

    ComputeStrain(elmtinfo,elmtsoln,_Strain);
    ComputeStressAndJacobian(InputParams,_Strain,_Stress,_Jac);

    _I.SetToIdentity();
    _devStress=_Stress-_I*(_Stress.Trace()/3.0);
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jac;

    Mate.Rank2Materials("F")=_F;
    Mate.Rank2Materials("PK1")=_Stress;//1-st Piola-Kirchhoff stress
    Mate.Rank2Materials("PK2")=_pk2;   //2-nd Piola-Kirchhoff stress

}
