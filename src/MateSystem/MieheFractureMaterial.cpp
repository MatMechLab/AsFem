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


void MieheFractureMaterial::InitMaterialProperties(const int &nDim,const Vector3d &gpCoord, const vector<double> &InputParams,
                                                   const vector<double> &gpU, const vector<double> &gpUdot,
                                                   const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUdot,
                                                   Materials &Mate) {
    // get rid of unused warning
    if(nDim||gpCoord(1)||InputParams.size()||gpU[0]||gpUdot[0]||
       gpGradU[0](1)||gpGradUdot[0](1)||Mate.ScalarMaterials.size()){
    }

    Mate.ScalarMaterials["Hist"]=0.0;
}
//***********************************************************
double MieheFractureMaterial::DegradationFun(const double &x) {
    return (1-x)*(1-x);
}
double MieheFractureMaterial::DegradationFunDeriv(const double &x) {
    return 2*(x-1);
}
//************************************************************
void MieheFractureMaterial::ComputeStrain(const int &nDim, const vector<Vector3d> &GradDisp, RankTwoTensor &strain) {
    // here we assume the first dof is d, then the following dofs are ux, uy and uz
    if(nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model only works for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(nDim==2){
        // DoFs: d ux uy
        GradU.SetFromGradU(GradDisp[2],GradDisp[3]);
    }
    else if(nDim==3){
        // DoFs: d ux uy uz
        GradU.SetFromGradU(GradDisp[2],GradDisp[3],GradDisp[4]);
    }
    strain=(GradU+GradU.Transpose())*0.5;
}
//************************************************************
void MieheFractureMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const RankTwoTensor &strain,const double &damage,
                                                    const Materials &MateOld, Materials &Mate) {
    double d;
    double g,dg;// for the degradation function
    const double k=1.0e-5; // for stabilizer

    const double lambda=InputParams[0];
    const double mu=InputParams[1];
    const double Gc=InputParams[2];
    const double L=InputParams[3];
    const double viscosity=InputParams[4];

    Mate.ScalarMaterials["viscosity"]=viscosity;
    Mate.ScalarMaterials["Gc"]=Gc;
    Mate.ScalarMaterials["L"]=L;

    int UseHist=0;
    if(InputParams.size()>=6){
        UseHist=static_cast<int>(InputParams[6-1]);
        if(UseHist<0) UseHist=0;
    }

    d=damage;
    g= DegradationFun(d);
    dg= DegradationFunDeriv(d);   // derivative of g

    I4Sym.SetToIdentitySymmetric4();

    Mate.Rank2Materials["strain"]=strain;

    ProjPos=strain.GetPostiveProjTensor();
    I4Sym.SetToIdentitySymmetric4();
    ProjNeg=I4Sym-ProjPos;

    // for the positive and negative strain
    EpsPos=ProjPos.DoubleDot(strain);
    EpsNeg=strain-EpsPos;

    double trEps,signpos,signneg;
    double psi,psipos,psineg;

    trEps=strain.Trace();
    psipos=0.5*lambda*BracketPos(trEps)*BracketPos(trEps)+mu*(EpsPos*EpsPos).Trace();
    psineg=0.5*lambda*BracketNeg(trEps)*BracketNeg(trEps)+mu*(EpsNeg*EpsNeg).Trace();
    psi=(g+k)*psipos+psineg;

    Mate.ScalarMaterials["Psi"]=psi;
    Mate.ScalarMaterials["PsiPos"]=psipos;
    Mate.ScalarMaterials["PsiNeg"]=psineg;


    I.SetToIdentity();
    StressPos=I*lambda*BracketPos(trEps)+EpsPos*2*mu;
    StressNeg=I*lambda*BracketNeg(trEps)+EpsNeg*2*mu;

    Mate.Rank2Materials["stress"]=StressPos*(g+k)+StressNeg;
    Mate.Rank2Materials["dstressdD"]=StressPos*dg;

    // for vonMises stress
    double trace;
    I.SetToIdentity();
    trace=Mate.Rank2Materials["stress"].Trace();
    DevStress=Mate.Rank2Materials["stress"]-I*(trace/3.0);
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    Mate.ScalarMaterials["vonMises"]=sqrt(1.5*DevStress.DoubleDot(DevStress));

    if(psipos>MateOld.ScalarMaterials.at("Hist")){
        Mate.ScalarMaterials["Hist"]=psipos;
        Mate.Rank2Materials["dHdstrain"]=StressPos;
    }
    else{
        Mate.ScalarMaterials["Hist"]=MateOld.ScalarMaterials.at("Hist");
        Mate.Rank2Materials["dHdstrain"].SetToZeros();
    }

    if(UseHist){
        Mate.ScalarMaterials["Hist"]=MateOld.ScalarMaterials.at("Hist");
        Mate.Rank2Materials["dHdstrain"].SetToZeros();
    }

    signpos=0.0;
    if(BracketPos(trEps)>0) signpos=1.0;

    signneg=0.0;
    if(BracketNeg(trEps)<0) signneg=1.0;

    Mate.Rank4Materials["jacobian"].SetToZeros();
    Mate.Rank4Materials["jacobian"]=(I.CrossDot(I)*lambda*signpos+ProjPos*2*mu)*(g+k)
            +I.CrossDot(I)*lambda*signneg+ProjNeg*2*mu;
}
//************************************************************
void MieheFractureMaterial::ComputeMaterialProperties(const double &t, const double &dt, const int &nDim,
                                                      const Vector3d &gpCoord, const vector<double> &InputParams,
                                                      const vector<double> &gpU, const vector<double> &gpUOld,
                                                      const vector<double> &gpUdot, const vector<double> &gpUdotOld,
                                                      const vector<Vector3d> &gpGradU,
                                                      const vector<Vector3d> &gpGradUOld,
                                                      const vector<Vector3d> &gpGradUdot,
                                                      const vector<Vector3d> &gpGradUdotOld, const Materials &MateOld,
                                                      Materials &Mate) {
    if(t||dt||nDim||gpCoord(1)||InputParams.size()||
       gpU[0]||gpUOld[0]||gpUdot[0]||gpUdotOld[0]||
       gpGradU[0](1)||gpGradUOld[0](1)||gpGradUdot[0](1)||gpGradUdotOld[0](1)||
       Mate.ScalarMaterials.size()||MateOld.ScalarMaterials.size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }

    // 1st dof: damage
    // 2nd dof: ux
    // 3rd dof: uy
    // 4th dof: uz

    ComputeStrain(nDim,gpGradU,Strain);
    ComputeConstitutiveLaws(InputParams,Strain,gpU[1],MateOld,Mate);
}