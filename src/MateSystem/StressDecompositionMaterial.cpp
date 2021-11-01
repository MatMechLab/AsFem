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
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model with stress decomposition
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/StressDecompositionMaterial.h"

void StressDecompositionMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    // get rid of unused warnings
    // **************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}
    Mate.Rank2Materials("stress").SetToZeros();
    Mate.ScalarMaterials("Hist")=0.0;
}
//***********************************************************
double StressDecompositionMaterial::DegradationFun(const double &x) {
    return (1-x)*(1-x);
}
double StressDecompositionMaterial::DegradationFunDeriv(const double &x) {
    return 2*(x-1);
}
//************************************************************
void StressDecompositionMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain){

    // here we assume the first dof is d, then the following dofs are ux, uy and uz
    if(elmtinfo.nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model only works for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(elmtinfo.nDim==2){
        // DoFs: d ux uy
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==3){
        // DoFs: d ux uy uz
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    Strain=(_GradU+_GradU.Transpose())*0.5;
}
//************************************************************
void StressDecompositionMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const double &damage,const RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian,const Materials &MateOld, Materials &Mate){
    
    double d;
    double g,dg;// for the degradation function
    const double k=1.0e-5; // for stabilizer

    double E=InputParams[1-1];
    double nu=InputParams[2-1];
    double Gc=InputParams[3-1];
    double L=InputParams[4-1];
    double viscosity=InputParams[5-1];
    double theta1,theta2,theta3;
    if(InputParams.size()<=5){
        theta1=0.0;theta2=0.0;theta3=0.0;
    }
    else{
        theta1=InputParams[6-1];
        theta2=InputParams[7-1];
        theta3=InputParams[8-1];
    }
    // here the theta must be degree, not rad
    _Rot.SetRotationTensorFromEulerAngle(theta1,theta2,theta3);
    _Rot.SetToIdentity();

    Mate.ScalarMaterials("viscosity")=viscosity;
    Mate.ScalarMaterials("Gc")=Gc;
    Mate.ScalarMaterials("L")=L;

    int UseHist=0;
    if(InputParams.size()>=9){
        UseHist=static_cast<int>(InputParams[9-1]);
        if(UseHist<0) UseHist=0;
    }

    d=damage;
    g= DegradationFun(d);
    dg= DegradationFunDeriv(d);   // derivative of g

    // fill the Cijkl tensor
    _Cijkl0.SetFromEandNu(E,nu);
    //_Cijkl=_Cijkl0.Rotate(_Rot);
    _Cijkl.SetFromEandNu(E,nu);

    // this stress is the temporary one!!!
    Stress=_Cijkl.DoubleDot(Strain);

    _ProjPos=Stress.GetPostiveProjTensor();
    _I4Sym.SetToIdentitySymmetric4();
    _ProjNeg=_I4Sym-_ProjNeg;

    // now we can get positive and negative stress
    _StressPos=_ProjPos.DoubleDot(Stress);
    _StressNeg=Stress-_StressPos;
    Stress=_StressPos*(g+k)+_StressNeg;
    // for the positive and negative elastic free energy
    _psipos=0.5*_StressPos.DoubleDot(Strain);
    _psineg=0.5*_StressNeg.DoubleDot(Strain);
    _psi=(g+k)*_psipos+_psineg;

    Mate.ScalarMaterials("Psi")=_psi;
    Mate.ScalarMaterials("PsiPos")=_psipos;
    Mate.ScalarMaterials("PsiNeg")=_psineg;

    Mate.Rank2Materials("dstressdD")=_StressPos*dg;

    Jacobian=_ProjPos.DoubleDot(_Cijkl)*(g+k)+_ProjNeg.DoubleDot(_Cijkl);
    Stress=Jacobian.DoubleDot(Strain);
    Mate.Rank2Materials("dstressdD")=(_ProjPos.DoubleDot(_Cijkl)).DoubleDot(Strain)*dg;

    Mate.ScalarMaterials("H")=0.0;// we use H instead of Hist for our element
                                  // in such a way, users can use the stagger solution!
    if(_psipos>MateOld.ScalarMaterials("Hist")){
        Mate.ScalarMaterials("Hist")=_psipos;
        Mate.Rank2Materials("dHdstrain")=_StressPos;
    }
    else{
        Mate.ScalarMaterials("Hist")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }

    Mate.ScalarMaterials("H")=Mate.ScalarMaterials("Hist");

    if(UseHist){
        Mate.ScalarMaterials("H")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }

}
//************************************************************
void StressDecompositionMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //***************************************************************
    // get rid of unused warnings 
    // **************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }
    else if(InputParams.size()>5 && InputParams.size()<8){
        MessagePrinter::PrintErrorTxt("if you want to give the euler angles, then please add them next to the viscosity coefficient, then, you need 8 parameters in total");
        MessagePrinter::PrintErrorTxt("the complete parameters list is: lambda, mu, Gc, L, viscosity, theta1, theta2, theta3, usehist-flag");
        MessagePrinter::AsFem_Exit();
    }

    // 1st dof: damage
    // 2nd dof: ux
    // 3rd dof: uy
    // 4th dof: uz

    ComputeStrain(elmtinfo,elmtsoln,_Strain);

    ComputeConstitutiveLaws(InputParams,elmtsoln.gpU[1],_Strain,_Stress,_Jacobian,MateOld,Mate);
    // for vonMises stress
    _I.SetToIdentity();
    _DevStress=_Stress-_I*(_Stress.Trace()/3.0);
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_DevStress.DoubleDot(_DevStress));

    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jacobian;

}
