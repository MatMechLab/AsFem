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
//+++ Date   : 2020.12.30
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"
#include "Utils/MathFuns.h"

void BulkMateSystem::MieheFractureMaterial(const int &nDim, const double &t, const double &dt,
                                          const vector<double> &InputParams, const Vector3d &gpCoord,
                                          const vector<double> &gpU, const vector<double> &gpV,
                                          const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradV,
                                          vector<double> &gpHist, const vector<double> &gpHistOld) {

    //*****************************************************************************
    //*** just to get rid of warnings, normal users dont need to do this
    //*****************************************************************************
    if(nDim||t||dt||InputParams.size()||
       gpCoord(1)||gpU.size()||gpV.size()||gpGradU.size()||gpGradV.size()||
       gpHist.size()||gpHistOld.size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }

    int UseHist=0;
    double d;
    double g,dg;// for the degradation function
    const double k=1.0e-6; // for stabilizer

    const double lambda=InputParams[0];
    const double mu=InputParams[1];
    const double Gc=InputParams[2];
    const double L=InputParams[3];
    const double viscosity=InputParams[4];


    _ScalarMaterials["viscosity"]=viscosity;
    _ScalarMaterials["Gc"]=Gc;
    _ScalarMaterials["L"]=L;


    UseHist=0;
    if(InputParams.size()>=6){
        UseHist=static_cast<int>(InputParams[6-1]);
        if(UseHist<0) UseHist=0;
    }


    d=gpU[1];
    g=(1-d)*(1-d);// degradation
    dg=2*(d-1);   // derivative of g

    RankTwoTensor Eps,EpsPos,EpsNeg,GradU;
    RankFourTensor ProjPos,ProjNeg;
    RankFourTensor I4Sym(RankFourTensor::InitIdentitySymmetric4);

    GradU.SetToZeros();
    if(nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model works only for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(nDim==2){
        GradU.SetFromGradU(gpGradU[2],gpGradU[3]);
    }
    else if(nDim==3){
        GradU.SetFromGradU(gpGradU[2],gpGradU[3],gpGradU[4]);
    }
    Eps=(GradU+GradU.Transpose())*0.5;

    _Rank2Materials["strain"]=Eps;

    ProjPos=Eps.GetPostiveProjTensor();
    I4Sym.SetToIdentitySymmetric4();
    ProjNeg=I4Sym-ProjPos;

    // for the positive and negative strain
    EpsPos=ProjPos.DoubleDot(Eps);
    EpsNeg=Eps-EpsPos;

    double trEps,signpos,signneg;
    double psi,psipos,psineg;

    trEps=Eps.Trace();
    psipos=0.5*lambda*BracketPos(trEps)*BracketPos(trEps)+mu*(EpsPos*EpsPos).Trace();
    psineg=0.5*lambda*BracketNeg(trEps)*BracketNeg(trEps)+mu*(EpsNeg*EpsNeg).Trace();
    psi=(g+k)*psipos+psineg;

    _ScalarMaterials["Psi"]=psi;
    _ScalarMaterials["PsiPos"]=psipos;
    _ScalarMaterials["PsiNeg"]=psineg;

    RankTwoTensor StressPos,StressNeg,I;

    I.SetToIdentity();
    StressPos=I*lambda*BracketPos(trEps)+EpsPos*2*mu;
    StressNeg=I*lambda*BracketNeg(trEps)+EpsNeg*2*mu;

    _Rank2Materials["stress"]=StressPos*(g+k)+StressNeg;
    _Rank2Materials["dstressdD"]=StressPos*dg;

    // for vonMises stress
    RankTwoTensor devStress;
    double trace;
    I.SetToIdentity();
    trace=_Rank2Materials["stress"].Trace();
    devStress=_Rank2Materials["stress"]-I*(trace/3.0);
    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    _ScalarMaterials["vonMises"]=sqrt(1.5*devStress.DoubleDot(devStress));

    if(psipos>gpHistOld[0]){
        _ScalarMaterials["Hist"]=psipos;
        _Rank2Materials["dHdstrain"]=StressPos;
    }
    else{
        _ScalarMaterials["Hist"]=gpHistOld[0];
        _Rank2Materials["dHdstrain"].SetToZeros();
    }

    if(UseHist){
        _ScalarMaterials["Hist"]=gpHistOld[0];
        _Rank2Materials["dHdstrain"].SetToZeros();
    }

    signpos=0.0;
    if(BracketPos(trEps)>0) signpos=1.0;

    signneg=0.0;
    if(BracketNeg(trEps)<0) signneg=1.0;

    _Rank4Materials["elasticity_tensor"].SetToZeros();
    _Rank4Materials["elasticity_tensor"]=(I.CrossDot(I)*lambda*signpos+ProjPos*2*mu)*(g+k)
                                         +I.CrossDot(I)*lambda*signneg+ProjNeg*2*mu;


}

