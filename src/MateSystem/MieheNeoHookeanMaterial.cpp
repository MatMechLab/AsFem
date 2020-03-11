//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::MieheNeoHookeanMaterial(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<5){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for phasefield fracture, 5 parameters are required   !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        E,nu,Gc,L,viscosity are expected for Miehe's model   !!!   ***\n");
        Msg_AsFem_Exit();
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
    const double EE=InputParams[0]; // Young's modulus
    const double nu=InputParams[1]; // Poisson ratio
    const double K=EE/(3*(1-2*nu)); // Bulk modulus
    const double G=EE/(2*(1+nu));

    _ScalarMaterials[0]=InputParams[4];// viscosity
    _ScalarMaterials[1]=InputParams[2];// Gc
    _ScalarMaterials[2]=InputParams[3];// L


    RankTwoTensor Stress(0.0),StressPos(0.0),StressNeg(0.0);
    RankTwoTensor F(0.0),Ft(0.0),C(0.0),Cinv(0.0),I(0.0),E(0.0);
    double J,J23,I1,I1bar;

    if(nDim==2){
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    I.SetToIdentity();
    F=_Rank2Materials[2]+I;
    J=F.Det();
    J23=pow(J,-2.0/3.0);//J23=J^{-2/3}
    
    Ft=F.Transpose();
    C=Ft*F;
    E=(C-I)*0.5;
    I1=C.Trace();
    I1bar=J23*I1;

    Cinv=C.Inverse();// C^{-1}

    // our total strain
    _Rank2Materials[0]=E;
    double d=gpU[nDim];
    if(d<0.0) d=0.0;
    if(d>1.0-1.0e-3) d=1.0;
    const double k=1.0e-5;// to avoid the zero stiffness matrix

    // for the fracture free energy
    double Psi,PsiPos,PsiNeg;

    if(J>=1.0){
        PsiPos=0.5*K*(J-1)*(J-1)+0.5*G*(I1bar-3.0);
        PsiNeg=0.0;
        Psi=(1-d)*(1-d)*PsiPos;

        StressPos=Cinv*K*(J-1)*J-Cinv*(G/3.0)*J23*I1+I*G*J23;
        StressNeg.SetToZeros();
        Stress=StressPos*((1-d)*(1-d)+k);
        // our constitutive matrix (Jacobian)
        _Rank4Materials[0]=(
            Cinv.CrossDot(Cinv)*K*(2*J-1)*J
           -Cinv.ODot(Cinv)*2*K*(J-1)*J 
           +(I.CrossDot(Cinv)*(-1.0/3.0)
            +Cinv.CrossDot(Cinv)*(I1/9.0)
            -Cinv.CrossDot(I)*(1.0/3.0)
            +Cinv.ODot(Cinv)*(I1/3.0))*2.0*G*J23
        )*((1-d)*(1-d)+k);
    }
    else{
        PsiPos=0.5*G*(I1bar-3.0);
        PsiNeg=0.5*K*(J-1)*(J-1);
        Psi=(1-d)*(1-d)*PsiPos+PsiNeg;
        StressPos=Cinv*(-G/3.0)*J23*I1+I*G*J23;
        StressNeg=Cinv*K*(J-1)*J;
        Stress=StressPos*((1-d)*(1-d)+k)+StressNeg;
        // our constitutive matrix (Jacobian)
        _Rank4Materials[0]=(
            I.CrossDot(Cinv)*(-1.0/3.0)
           +Cinv.CrossDot(Cinv)*(I1/9.0)
           -Cinv.CrossDot(I)*(1.0/3.0)
           +Cinv.ODot(Cinv)*(I1/3.0)
            )*((1-d)*(1-d)+k)*2.0*G*J23
        +Cinv.CrossDot(Cinv)*K*(2*J-1)*J
        -Cinv.ODot(Cinv)*2*K*(J-1)*J;
    }

    // out final stress
    _Rank2Materials[1]=Stress;
    _Rank2Materials[2]=StressPos*(-2+2*d);// dstress/dD




    // calculate H, and update the history variable
    _Rank2Materials[3].SetToZeros();// store the dH/dstrain rank2 tensor
    if(PsiPos>gpHistOld[0]){
        gpHist[0]=PsiPos;
        _Rank2Materials[3]=StressPos;
    }
    else{
        gpHist[0]=gpHistOld[0];
        _Rank2Materials[3].SetToZeros();
    }
    _Rank2Materials[3].SetToZeros();


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