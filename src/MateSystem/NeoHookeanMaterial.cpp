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
//+++ Date   : 2021.04.11
//+++ Purpose: Implement the calculation of neo-hookean type
//+++          hyperlelastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/NeoHookeanMaterial.h"


void NeoHookeanMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}
}

//*********************************************************
void NeoHookeanMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, RankTwoTensor &Strain){
    // here we calculate the deformation gradient tensor as well as euler-lagrange strain tensor
    if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    _I.SetToIdentity();
    _F=_GradU+_I;// F=I+U_{i,j}
    _C=_F.Transpose()*_F; //C=F^T*F
    _Cinv=_C.Inverse();
    Strain=(_C-_I)*0.5;

}
//****************************************************************
void NeoHookeanMaterial::ComputeStressAndJacobian(const vector<double> &InputParams, const RankTwoTensor &Strain, RankTwoTensor &Stress, RankFourTensor &Jacobian){
    // here Strain is E
    if(Strain(1,1)){}// get rid of unused warning

    
    double EE=InputParams[0];
    double nu=InputParams[1];

    double lambda=EE*nu/((1+nu)*(1-2*nu));
    double mu=EE/(2*(1+nu));
    double J=_F.Det();

    Stress=(_I-_Cinv)*mu+_Cinv*lambda*log(J);
    Jacobian=_Cinv.ODot(_Cinv)*(mu-lambda*log(J))*2
            +_Cinv.CrossDot(_Cinv)*lambda;

}
//***************************************************************
void NeoHookeanMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    if(InputParams.size()<2){
        MessagePrinter::PrintErrorTxt("for the NeoHookeanMaterial, two parameters are required, you need to give: E and nu");
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

}
