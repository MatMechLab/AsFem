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
//+++ Date   : 2021.04.09
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/MieheFractureMaterial.h"

void MieheFractureMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    // get rid of unused warnings
    // **************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.ScalarMaterials("Hist")=0.0;
    Mate.Rank2Materials("stress").SetToZeros();
}
//***********************************************************
double MieheFractureMaterial::DegradationFun(const double &x) {
    return (1-x)*(1-x);
}
double MieheFractureMaterial::DegradationFunDeriv(const double &x) {
    return 2*(x-1);
}
//************************************************************
void MieheFractureMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain){

    // here we assume the first dof is d, then the following dofs are ux, uy and uz
    if(elmtinfo.nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model only works for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(elmtinfo.nDim==2){
        // DoFs: d ux uy
        GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==3){
        // DoFs: d ux uy uz
        GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    Strain=(GradU+GradU.Transpose())*0.5;
}
//************************************************************
void MieheFractureMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const double &damage,const RankTwoTensor &strain,RankTwoTensor &Stress,RankFourTensor &Jacobian,const Materials &MateOld, Materials &Mate){
    
    double d;
    double g,dg;// for the degradation function
    const double k=5.0e-5; // for stabilizer

    const double lambda=InputParams[0];
    const double mu=InputParams[1];
    const double Gc=InputParams[2];
    const double L=InputParams[3];
    const double viscosity=InputParams[4];

    Mate.ScalarMaterials("viscosity")=viscosity;
    Mate.ScalarMaterials("Gc")=Gc;
    Mate.ScalarMaterials("L")=L;

    int UseHist=0;
    if(InputParams.size()>=6){
        UseHist=static_cast<int>(InputParams[6-1]);
        if(UseHist<0) UseHist=0;
    }

    d=damage;
    g= DegradationFun(d);
    dg= DegradationFunDeriv(d);   // derivative of g

    I4Sym.SetToIdentitySymmetric4();


    ProjPos=Strain.GetPositiveProjTensor();
    I4Sym.SetToIdentitySymmetric4();
    ProjNeg=I4Sym-ProjPos;

    // for the positive and negative strain
    EpsPos=ProjPos.DoubleDot(strain);
    EpsNeg=Strain-EpsPos;

    double trEps,signpos,signneg;
    double psi,psipos,psineg;

    trEps=Strain.Trace();
    psipos=0.5*lambda*BracketPos(trEps)*BracketPos(trEps)+mu*(EpsPos*EpsPos).Trace();
    psineg=0.5*lambda*BracketNeg(trEps)*BracketNeg(trEps)+mu*(EpsNeg*EpsNeg).Trace();
    psi=(g+k)*psipos+psineg;

    Mate.ScalarMaterials("Psi")=psi;
    Mate.ScalarMaterials("PsiPos")=psipos;
    Mate.ScalarMaterials("PsiNeg")=psineg;


    I.SetToIdentity();
    StressPos=I*lambda*BracketPos(trEps)+EpsPos*2*mu;
    StressNeg=I*lambda*BracketNeg(trEps)+EpsNeg*2*mu;

    Stress=StressPos*(g+k)+StressNeg;
    Mate.Rank2Materials("dstressdD")=StressPos*dg;

    Mate.ScalarMaterials("H")=0.0;// we use H instead of Hist for our element
                                  // in such a way, users can use the stagger solution!
    if(psipos>MateOld.ScalarMaterials("Hist")){
        Mate.ScalarMaterials("Hist")=psipos;
        Mate.Rank2Materials("dHdstrain")=StressPos;
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

    signpos=0.0;
    if(BracketPos(trEps)>0) signpos=1.0;

    signneg=0.0;
    if(BracketNeg(trEps)<0) signneg=1.0;

    Jacobian.SetToZeros();
    Jacobian=(I.OTimes(I)*lambda*signpos+ProjPos*2*mu)*(g+k)
        +I.OTimes(I)*lambda*signneg+ProjNeg*2*mu;

}
//************************************************************
void MieheFractureMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //***************************************************************
    // get rid of unused warnings 
    // **************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }

    // 1st dof: damage
    // 2nd dof: ux
    // 3rd dof: uy
    // 4th dof: uz

    ComputeStrain(elmtinfo,elmtsoln,Strain);


    ComputeConstitutiveLaws(InputParams,elmtsoln.gpU[1],Strain,Stress,Jacobian,MateOld,Mate);
    // for vonMises stress
    double trace;
    I.SetToIdentity();
    trace=Stress.Trace();
    DevStress=Stress-I*(trace/3.0);
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*DevStress.DoubleDot(DevStress));

    Mate.Rank2Materials("strain")=Strain;
    Mate.Rank2Materials("stress")=Stress;
    Mate.Rank4Materials("jacobian")=Jacobian;

    // used in AllenCahn fracture element
    Mate.ScalarMaterials("dFdD")=elmtsoln.gpU[1];
    Mate.ScalarMaterials("d2FdD2")=1.0;
    Mate.ScalarMaterials("M")=1.0/InputParams[5-1];// M=1.0/viscosity

}
