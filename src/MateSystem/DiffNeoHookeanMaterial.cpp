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
//+++ Date   : 2021.11.14
//+++ Purpose: Implement the calculation of diffusion coupled
//+++           neo-hookean type hyperlelastic material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/DiffNeoHookeanMaterial.h"


void DiffNeoHookeanMaterial::InitMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();
    Mate.Rank2Materials("F").SetToIdentity();
    Mate.Rank2Materials("PK1").SetToZeros();
    Mate.Rank2Materials("PK2").SetToZeros();
}

//*********************************************************
void DiffNeoHookeanMaterial::ComputeDeformationGradientTensor(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &F){
    // here we calculate the deformation gradient tensor as well as euler-lagrange strain tensor
    // DoFs:
    //   1) concentration
    //   2) ux
    //   3) uy
    //   4) uz
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[2]);
    }
    else if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    _I.SetToIdentity();
    F=_GradU+_I;// F=I+U_{i,j}
}
//****************************************************************
void DiffNeoHookeanMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const LocalElmtSolution &soln,const RankTwoTensor &F,RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian){
    double E=InputParams[3-1];
    double nu=InputParams[4-1];
    double PsiE,Omega,K,G;
    double Js,Je,Je23,dJsdC,C;
    double I1,I1bar;

    K=E/(3*(1-2*nu));// bulk modulus
    G=E/(2*(1+nu));  // shear modulus

    Omega=InputParams[2-1];
    C=soln.gpU[1];

    // the elastic free energy is:
    // psi=Js*psie
    // for the swelling part from diffusion
    Js=1.0+Omega*C;
    dJsdC=Omega;
    _I.SetToIdentity();
    _Fc=_I*cbrt(Js);

    // for the elastic deformation contribution
    _Fe=F*_Fc.Inverse();
    Je=_Fe.Det();
    Je23=pow(Je,-2.0/3.0);
    _Ce=_Fe.Transpose()*_Fe;
    _CeInv=_Ce.Inverse();
    I1=_Ce.Trace();
    I1bar=I1*Je23;

    // now we can calculate the stress and jacobian
    Strain=(_Ce-_I)*0.5;

    Stress=_CeInv*0.5*K*(Je*Je-1)
        -_CeInv*(G/3.0)*I1bar
        +_I*G*Je23;

    Jacobian.SetToZeros();
    Jacobian=_CeInv.OTimes(_CeInv)*K*Je*Je-_CeInv.ODot(_CeInv)*K*(Je*Je-1)
        +(
         _CeInv.OTimes(_CeInv)*(I1/3.0)
        -_CeInv.OTimes(_I)
        -_I.OTimes(_CeInv)
        +_CeInv.ODot(_CeInv)*I1
         )*(2*G/3.0)*Je23;

    // elastic free energy density is:
    //  PsiE=0.5*K(0.5*(Je^2-1)-log(Je))+0.5*G*(I1bar-3.0)
    PsiE=0.5*K*(0.5*(Je*Je-1)-log(Je))+0.5*G*(I1bar-3.0);
    _Mu=dJsdC*PsiE;
    _dStressdC.SetToZeros();
    _dMudStrain.SetToZeros();
    _dMudStrain=Stress*dJsdC;

}
//***************************************************************
void DiffNeoHookeanMaterial::ComputeMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,
                                           const Materials &MateOld,Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    if(InputParams.size()<4){
        MessagePrinter::PrintErrorTxt("for the DiffNeoHookeanMaterial, four parameters are required, you need to give: D, Omega, E, and nu");
        MessagePrinter::AsFem_Exit();
    }
    Mate.ScalarMaterials("D")=InputParams[1-1];
    Mate.ScalarMaterials("Omega")=InputParams[2-1];

    ComputeDeformationGradientTensor(elmtinfo,elmtsoln,_F);
    ComputeConstitutiveLaws(InputParams,elmtsoln,_F,_Strain,_Stress,_Jac);

    Mate.Rank2Materials("F")=_F;
    Mate.Rank2Materials("PK1")=_F*_Stress;//1-st Piola-Kirchhoff stress 
    Mate.Rank2Materials("PK2")=_Stress;   //2-nd Piola-Kirchhoff stress

    _Stress=_Fe*_Stress;// we use 1-st Piola-Kirchhoff stress for the calculation

    _devStress=_Stress.Dev();
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jac;

    Mate.VectorMaterials("GradSigmaH")=0.0;
    Mate.ScalarMaterials("dSigmaHdC")=0.0;
    

}
