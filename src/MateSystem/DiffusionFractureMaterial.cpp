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
//+++ Date   : 2021.11.14
//+++ Purpose: Implement the calculation of diffusion coupled
//+++          phase-field fracture material
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/DiffusionFractureMaterial.h"
#include "Utils/MathFuns.h"

void DiffusionFractureMaterial::InitMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();
    Mate.ScalarMaterials("Hist")=0.0;
}

//*********************************************************
void DiffusionFractureMaterial::ComputeDeformationGradientTensor(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &F){
    // here we calculate the deformation gradient tensor as well as euler-lagrange strain tensor
    // DoFs:
    //   1) c 
    //   2) d
    //   2) ux
    //   3) uy
    //   4) uz
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3],elmtsoln.gpGradU[4],elmtsoln.gpGradU[5]);
    }
    F=_GradU;// for small strain case, we set F as the disp gradient, 
             //it is not the deformation gradient tensor !!!
}
//****************************************************************
void DiffusionFractureMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const LocalElmtSolution &soln,const RankTwoTensor &F,RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian){
    double E=InputParams[6-1];
    double nu=InputParams[7-1];
    double Omega,lambda,mu;
    double C,d;
    double g,dg;

    lambda=E*nu/((1+nu)*(1-2*nu));// lame 1st constant
    mu=E/(2*(1+nu));              // shear modulus

    Omega=InputParams[2-1];
    C=soln.gpU[1];
    d=soln.gpU[2];

    _I.SetToIdentity();
    // for our total small strain
    Strain=(F+F.Transpose())*0.5;// F is set as GradU in 'ComputeDeformationGradientTensor'
    
    // for the diffusion induced eigen strain
    _EigenStrain=_I*Omega*C/3.0;
    _dEigenStraindC=_I*Omega/3.0;

    // for elastic strain
    _MechStrain=Strain-_EigenStrain;
    _dMechStraindC=_dEigenStraindC*-1.0;

    // for degradation function and its derivative
    g= DegradationFun(d);
    dg= DegradationFunDeriv(d);   // derivative of g

    _ProjPos=_MechStrain.GetPositiveProjTensor();
    _I4Sym.SetToIdentitySymmetric4();
    _ProjNeg=_I4Sym-_ProjPos;

    // for the positive and negative strain
    _EpsPos=_ProjPos.DoubleDot(_MechStrain);
    _EpsNeg=_MechStrain-_EpsPos;

    _dEpsPosdC=_ProjPos.DoubleDot(_dMechStraindC);
    _dEpsNegdC=_dMechStraindC-_dEpsPosdC;

    double trEps,signpos,signneg,dtrEpsdC;

    trEps=_MechStrain.Trace();
    dtrEpsdC=_dMechStraindC.Trace();
    _PsiPos=0.5*lambda*BracketPos(trEps)*BracketPos(trEps)+mu*(_EpsPos*_EpsPos).Trace();
    _PsiNeg=0.5*lambda*BracketNeg(trEps)*BracketNeg(trEps)+mu*(_EpsNeg*_EpsNeg).Trace();

    _I.SetToIdentity();
    _StressPos=_I*lambda*BracketPos(trEps)+_EpsPos*2*mu;
    _StressNeg=_I*lambda*BracketNeg(trEps)+_EpsNeg*2*mu;

    _dStressPosdC=_I*lambda*BracketPos(dtrEpsdC)+_dEpsPosdC*2*mu;
    _dStressNegdC=_I*lambda*BracketNeg(dtrEpsdC)+_dEpsNegdC*2*mu;

    Stress=_StressPos*(g+1.0e-5)+_StressNeg;

    _dStressdC=_dStressPosdC*(g+1.0e-5)+_dStressNegdC;
    _dStressdD=_StressPos*dg;

    _dSigmaHdC=_dStressdC.Trace()/3.0;
    _GradSigmaH=soln.gpGradU[1]*_dSigmaHdC;

    signpos=0.0;
    if(BracketPos(trEps)>0) signpos=1.0;

    signneg=0.0;
    if(BracketNeg(trEps)<0) signneg=1.0;

    Jacobian.SetToZeros();
    Jacobian=(_I.OTimes(_I)*lambda*signpos+_ProjPos*2*mu)*(g+1.0e-5)
        +_I.OTimes(_I)*lambda*signneg+_ProjNeg*2*mu;

}
//***************************************************************
void DiffusionFractureMaterial::ComputeMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,
                                           const Materials &MateOld,Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    if(InputParams.size()<7){
        MessagePrinter::PrintErrorTxt("for the DiffusionFractureMaterial, four parameters are required, you need to give: D, Omega, Gc, L0, Viscosity, E, and nu");
        MessagePrinter::AsFem_Exit();
    }
    Mate.ScalarMaterials("D")=InputParams[1-1];
    Mate.ScalarMaterials("Omega")=InputParams[2-1];
    Mate.ScalarMaterials("Gc")=InputParams[3-1];
    Mate.ScalarMaterials("L")=InputParams[4-1];
    Mate.ScalarMaterials("viscosity")=InputParams[5-1];


    ComputeDeformationGradientTensor(elmtinfo,elmtsoln,_F);
    ComputeConstitutiveLaws(InputParams,elmtsoln,_F,_Strain,_Stress,_Jac);

    if(_PsiPos>MateOld.ScalarMaterials("Hist")){
        Mate.ScalarMaterials("Hist")=_PsiPos;
        Mate.Rank2Materials("dHdstrain")=_StressPos;
    }
    else{
        Mate.ScalarMaterials("Hist")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }

    if(InputParams.size()>7){
        UseHist=static_cast<int>(InputParams[8-1]);
        if(UseHist<0) UseHist=0;
    }

    Mate.ScalarMaterials("H")=Mate.ScalarMaterials("Hist");
    if(UseHist){
        Mate.ScalarMaterials("H")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }

    _devStress=_Stress.Dev();
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jac;

    Mate.Rank2Materials("dStressdC")=_dStressdC;
    Mate.Rank2Materials("dStressdD")=_dStressdD;

    Mate.VectorMaterials("GradSigmaH")=_GradSigmaH;
    Mate.ScalarMaterials("dSigmaHdC")=_dSigmaHdC;

    // used in AllenCahn fracture element
    Mate.ScalarMaterials("dFdD")=elmtsoln.gpU[1];
    Mate.ScalarMaterials("d2FdD2")=1.0;
    Mate.ScalarMaterials("M")=1.0/InputParams[5-1];// M=1.0/viscosity

}
