//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::UserMaterial2(const int &nDim,const double &t,const double &dt,
                        const vector<double> InputParams,
                        const Vector3d &gpCoord,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                        vector<double> &gpHist,const vector<double> &gpHistOld){
    if(nDim){}
    if(t||dt){}
    if(gpCoord(0)){}
    if(gpU[0]){}
    if(gpV[0]){}
    if(gpGradU[0](0)){}
    if(gpGradV[0](0)){}
    if(gpHist[0]){}
    if(gpHistOld[0]){}

    if(InputParams.size()<7){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for phasefield fracture, 7 parameters are required   !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        Emin,Emax,nu,Gc,L,viscosity are expected for Miehe's model!***\n");
        Msg_AsFem_Exit();
    }
    // we use the sixth one to indicate wether we use the history one(stager or not)
    int UseHist=0;
    if(InputParams.size()>=6){
        if(int(InputParams[6-1])==1){
            UseHist=1;
        }
        else if(int(InputParams[6-1])==0){
            UseHist=0;
        }
        UseHist=0;
    }

    //*********************************************
    //*** IMPORTANT!!!
    //*** in this model, d=0--->for undamaged case
    //***                d=1--->for damaged   case
    //**********************************************

    //*******************************************
    //*** material design
    //*** viscosity=ScalarMaterials[0];
    //*** Gc=ScalarMaterials[1];
    //*** L=ScalarMaterials[2];
    //*** H=Hist[0];
    //*******************************************
    const double Emin=InputParams[0];
    const double Emax=InputParams[1];
    double dE=(Emax-Emin)/2.0;
    double EE=(Emin+Emax)/2.0+dE*sin(2*gpCoord(1)*gpCoord(2));
    double nu=InputParams[2]; // Poisson ratio
    double lambda=EE*nu/((1+nu)*(1-2*nu));// lame const
    double mu=EE/(2*(1+nu));

    _ScalarMaterials[0]=InputParams[5];// viscosity
    _ScalarMaterials[1]=InputParams[3];// Gc
    _ScalarMaterials[2]=InputParams[4];// L


    RankTwoTensor Stress(0.0),StressPos(0.0),StressNeg(0.0);
    RankTwoTensor Strain(0.0),StrainPos(0.0),StrainNeg(0.0);
    
    if(nDim==2){
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    Strain=(_Rank2Materials[2]+_Rank2Materials[2].Transpose())*0.5;
    // our total strain
    _Rank2Materials[0]=Strain;
    
    RankTwoTensor eigvec(0.0);
    double eigval[3];

    RankFourTensor ProjPos=Strain.CalcPostiveProjTensor(eigval,eigvec);
    RankFourTensor I4Sym(RankFourTensor::InitIdentitySymmetric4);
    RankFourTensor ProjNeg=I4Sym-ProjPos;

    StrainPos=ProjPos.DoubleDot(Strain);
    StrainNeg=Strain-StrainPos;

    double StrainTrace=Strain.Trace();

    double TrPos= (abs(StrainTrace)+StrainTrace)*0.5;
    double TrNeg=-(abs(StrainTrace)-StrainTrace)*0.5;

    // now we can split the positive and negative stress
    RankTwoTensor I(0.0);
    I.SetToIdentity();// Unity tensor
    StressPos=I*lambda*TrPos+StrainPos*2.0*mu;
    StressNeg=I*lambda*TrNeg+StrainNeg*2.0*mu;
    // now we can have the final stress
    double d=gpU[nDim];
    if(d<1.0e-2) d=1.0e-2;
    if(d>1.0-1.0e-2) d=1.0-1.0e-2;
    const double k=5.0e-4; // to avoid the zero stiffness matrix
    Stress=((1-d)*(1-d)+k)*StressPos+StressNeg;
    // store the stress and dstress/dd in rank2 material
    _Rank2Materials[1]=Stress;
    _Rank2Materials[2]=(-2+2*d)*StressPos;//dStress/dD
    // for our constitutive law
    double SignPos,SignNeg;
    SignPos=0.0;
    if(StrainTrace>=0.0) SignPos=1.0;
    SignNeg=0.0;
    if(StrainTrace<=0.0) SignNeg=1.0;
    _Rank4Materials[0]=(I.CrossDot(I)*lambda*SignPos+ProjPos*2*mu)*((1-d)*(1-d)+k)
                      +(I.CrossDot(I)*lambda*SignNeg+ProjNeg*2*mu);


    // for the fracture free energy
    double Psi,PsiPos,PsiNeg;
    PsiPos=0.5*lambda*TrPos*TrPos+mu*((StrainPos*StrainPos).Trace());
    PsiNeg=0.5*lambda*TrNeg*TrNeg+mu*((StrainNeg*StrainNeg).Trace());
    Psi=(1-d)*(1-d)*PsiPos+PsiNeg;


    // calculate H, and update the history variable
    _Rank2Materials[3].SetToZeros();// store the dH/dstrain rank2 tensor
    if(PsiPos>=gpHistOld[0]){
        gpHist[0]=PsiPos;
        _Rank2Materials[3]=StressPos;
    }
    else{
        gpHist[0]=gpHistOld[0];
        _Rank2Materials[3].SetToZeros();
    }
    
    //*************************************
    if(UseHist){
        // If we use the history state variable, then the stager solution is used
        _ScalarMaterials[20]=gpHistOld[0];
        _Rank2Materials[5].SetToZeros();
    }
    else{
        // Or we use the fully coupled case
        _ScalarMaterials[20]=gpHist[0];
        _Rank2Materials[5]=_Rank2Materials[3];
    }


    // use the fourth one to calculate the vonMises stress
    _Rank2Materials[4]=_Rank2Materials[1];
    // stress deviator tensor
    _Rank2Materials[4](1,1)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[4](2,2)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[4](3,3)-=(1.0/3.0)*_Rank2Materials[1].Trace();



    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _ScalarMaterials[3]=sqrt(1.5*_Rank2Materials[4].DoubleDot(_Rank2Materials[4]));
    // hydrostatic stress
    _ScalarMaterials[4]=_Rank2Materials[1].Trace()/3.0;
    if(nDim==2){
        _ScalarMaterials[ 5]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[ 6]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[ 7]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[ 8]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[ 9]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[10]=_Rank2Materials[0](1,2);//strxy

        _ScalarMaterials[11]=Psi;
        _ScalarMaterials[12]=PsiPos;
        _ScalarMaterials[13]=PsiNeg;
    }
    else if(nDim==3){
        _ScalarMaterials[ 5]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[ 6]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[ 7]=_Rank2Materials[1](3,3);//sigzz
        _ScalarMaterials[ 8]=_Rank2Materials[1](2,3);//sigyz
        _ScalarMaterials[ 9]=_Rank2Materials[1](1,3);//sigzx
        _ScalarMaterials[10]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[11]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[12]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[13]=_Rank2Materials[0](3,3);//strzz
        _ScalarMaterials[14]=_Rank2Materials[0](2,3);//stryz
        _ScalarMaterials[15]=_Rank2Materials[0](1,3);//strzx
        _ScalarMaterials[16]=_Rank2Materials[0](1,2);//strxy

        _ScalarMaterials[17]=Psi;
        _ScalarMaterials[18]=PsiPos;
        _ScalarMaterials[19]=PsiNeg;
    }

}