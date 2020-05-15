//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::NeoHookeanMaterial(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<3){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for neo-hookean mate, three parameters are required  !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        incompressive case needs E, nu, comp=0, pressure     !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        compressive case needs E,nu, comp=1                  !!!   ***\n");
        Msg_AsFem_Exit();
    }


    double EE=InputParams[0];
    double nu=InputParams[1];

    double lambda=EE*nu/((1+nu)*(1-2*nu));
    double mu=EE/(2*(1+nu));
    double Pressure=0.0;
    int IsCompressive=1;

    if(InputParams.size()>=3){
        if(static_cast<int>(InputParams[3-1])==1){
            IsCompressive=1;
        }
        else if(static_cast<int>(InputParams[3-1])==0){
            IsCompressive=0;
            Pressure=0.0;
            if(InputParams.size()<4){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for neo-hookean mate, 4 parameters are required      !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        incompressive case need: E, nu,comp=0,pressure       !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                Pressure=InputParams[4-1];
            }
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: for neo-hookean mate, three parameters are required  !!!   ***\n");
            PetscPrintf(PETSC_COMM_WORLD,"***        E and nu, comp are expected                          !!!   ***\n");
            PetscPrintf(PETSC_COMM_WORLD,"***        comp=1 for compressive, comp=0 for incompressive     !!!   ***\n");
            Msg_AsFem_Exit();
        }
    }
    //***********************************************
    //*** the related formulas are taken from:
    //*** https://en.wikipedia.org/wiki/Neo-Hookean_solid
    //************************************************
    //*** for model related parameters
    double W;
    double I1,I1bar;
    double I3;
    double J;
    RankTwoTensor F(0.0),Ft(0.0),C(0.0),Cinv(0.0);
    RankTwoTensor E(0.0),I(0.0);
    RankTwoTensor Stress(0.0),Strain(0.0);
    RankFourTensor Jac(0.0);


    //*******************************************
    //*** Firstly, we calculate the finite strain
    //*******************************************
    _Rank2Materials[2].SetToZeros();
    
    I.SetToIdentity();
    
    if(nDim==2){
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    
    // for the deformation gradient
    F=_Rank2Materials[2]+I;// F=I+U_{i,j}
    J=F.Det();
    // for the right Cauchy-Green tensor
    C=F.Transpose()*F;//C=F^tF
    Cinv=C.Inverse();
    I1=C.Trace();
    I3=C.Det();
    I1bar=I1/cbrt(I3);
    // for Green-Lagrange strain
    E=(C-I)*0.5;


    //******************************
    //*** our finite strain
    //******************************
    _Rank2Materials[0]=E;
    if(!IsCompressive){
        // for incompressive neo-hookean model
        // Pressure=3*lambda*(J-1);
        W=0.5*mu*(I1bar-3.0)+Pressure*(J-1);
        Stress=(I-Cinv*I1/3.0)*mu*(1.0/cbrt(I3))
              +Cinv*Pressure*J;
        Jac=( I.CrossDot(Cinv)*(-1.0/3.0)
             +Cinv.CrossDot(Cinv)*(I1/9.0)
             +Cinv.CrossDot(I)*(-1.0/3.0)
             +Cinv.ODot(Cinv)*(I1/3.0))*2*mu*(1.0/cbrt(I3))
           +(Cinv.CrossDot(Cinv)-Cinv.ODot(Cinv)*2)*Pressure*J;
    }
    else{
        // for compressive neo-hookean model
        W=0.5*mu*(I1-2-2*log(J))+0.5*lambda*(J-1)*(J-1);
        Stress=(I-Cinv)*mu+Cinv*lambda*(J-1)*J;
        Jac=Cinv.ODot(Cinv)*mu*2
           +Cinv.CrossDot(Cinv)*lambda*(2*J-1)*J
           -Cinv.ODot(Cinv)*lambda*(J-1)*J*2;
    }

    //******************************
    //*** our stress
    //******************************

    _Rank2Materials[1]=Stress;

    
    //******************************
    //*** now we calculate our jacobian matrix
    //******************************
    _Rank4Materials[0]=Jac;


    // use the fourth one to calculate the vonMises stress
    _Rank2Materials[3]=_Rank2Materials[1];
    // stress deviator tensor
    _Rank2Materials[3](1,1)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[3](2,2)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[3](3,3)-=(1.0/3.0)*_Rank2Materials[1].Trace();



    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _ScalarMaterials[0]=sqrt(1.5*_Rank2Materials[3].DoubleDot(_Rank2Materials[3]));
    // hydrostatic stress
    _ScalarMaterials[1]=_Rank2Materials[1].Trace()/3.0;
    if(nDim==2){
        _ScalarMaterials[2]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[3]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[4]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[5]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[6]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[7]=_Rank2Materials[0](1,2);//strxy
        _ScalarMaterials[8]=W;
    }
    else if(nDim==3){
        _ScalarMaterials[2]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[3]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[4]=_Rank2Materials[1](3,3);//sigzz
        _ScalarMaterials[5]=_Rank2Materials[1](2,3);//sigyz
        _ScalarMaterials[6]=_Rank2Materials[1](1,3);//sigzx
        _ScalarMaterials[7]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[ 8]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[ 9]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[10]=_Rank2Materials[0](3,3);//strzz
        _ScalarMaterials[11]=_Rank2Materials[0](2,3);//stryz
        _ScalarMaterials[12]=_Rank2Materials[0](1,3);//strzx
        _ScalarMaterials[13]=_Rank2Materials[0](1,2);//strxy
        _ScalarMaterials[14]=W;
    }

}