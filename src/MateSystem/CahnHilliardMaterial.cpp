//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::CahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for cahn hilliard mate, three parameters are required!!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        D, chi, kappa are expected                           !!!   ***\n");
        Msg_AsFem_Exit();
    }
    //****************************************
    //*** D, chi, kappa, free energy choice
    //*** for free-energy-choice=1(or only D, chi, kappa are given, then mate=1)
    //*** for free-energy-choice=2
    //***    except D,chi, kappa, mate, Calpha, Cbeta
    //****************************************
    double F=0,dFdc=0,d2Fdc2=0;
    double Source=0,dSourcedC=0,dSourcedMu=0;
    double Calpha=0,Cbeta=0;
    double D,M,dMdc;
    double Chi,Kappa;
    double ReacRate=0.0;
    int Mate=1;



    double c=gpU[0];// if you use cahn-hilliard,
                    // must make sure the first dof is concentration, not mu!!!
    double mu=gpU[1];
    double tol=1.0e-4;
    if(c<tol) c=tol;
    if(c>1.0-tol) c=1.0-tol;

    D=InputParams[1-1];
    Chi=InputParams[2-1];
    Kappa=InputParams[3-1];
    Mate=1;
    if(InputParams.size()>=4){
        Mate=static_cast<int>(InputParams[4-1]);
        if(Mate==2){
            if(InputParams.size()<6){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for cahn hilliard mate, three parameters are required!!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        D, chi, kappa, mate, calpha, cbeta are expected      !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                Calpha=InputParams[5-1];
                Cbeta=InputParams[6-1];
            }
        }
        else if(Mate==3){
            if(InputParams.size()<7){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for cahn hilliard mate, three parameters are required!!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        D, chi, kappa, mate, calpha, cbeta,K are expected    !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                Calpha=InputParams[5-1];
                Cbeta=InputParams[6-1];
                ReacRate=InputParams[7-1];
            }
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported material choice for CahnHilliard material!!!   ***\n");
            Msg_AsFem_Exit();
        }
    }

    M=D*c*(1-c);
    dMdc=D*(1-2*c);
    if(Mate==1){
        // for binary-mixed free energy and its derivative
        F=c*log(c)+(1-c)*log(1-c)+Chi*c*(1-c);
        dFdc=log(c)-log(1-c)+Chi*(1-2*c);
        d2Fdc2=1.0/c+1.0/(1.0-c)-Chi*2.0;
        Source=0.0;dSourcedC=0.0;dSourcedMu=0.0;
    }
    else if(Mate==2){
        // For the polynomial case
        F=(c-Calpha)*(c-Calpha)*(c-Cbeta)*(c-Cbeta);
        dFdc=2*(c-Calpha)*(c-Cbeta)*(2*c-Calpha-Cbeta);
        d2Fdc2=2*(Calpha*Calpha+4*Calpha*Cbeta-6*Calpha*c+Cbeta*Cbeta-6*Cbeta*c+6*c*c);
        Source=0.0;dSourcedC=0.0;dSourcedMu=0.0;
    }
    else if(Mate==3){
        F=(c-Calpha)*(c-Calpha)*(c-Cbeta)*(c-Cbeta);
        dFdc=2*(c-Calpha)*(c-Cbeta)*(2*c-Calpha-Cbeta);
        d2Fdc2=2*(Calpha*Calpha+4*Calpha*Cbeta-6*Calpha*c+Cbeta*Cbeta-6*Cbeta*c+6*c*c);
        Source=ReacRate*(1-c)*(exp(mu/2.0)-exp(-mu/2.0));
        dSourcedC=-ReacRate*(exp(mu/2.0)-exp(-mu/2.0));
        dSourcedMu=ReacRate*(1-c)*(exp(mu/2.0)+exp(-mu/2.0))*0.5;
    }


    _ScalarMaterials[0]=M;// M=D*c*(1-c)
    _ScalarMaterials[1]=dMdc;// dM/dc
    _ScalarMaterials[2]=Chi;// Chi
    _ScalarMaterials[3]=Kappa;// Kappa
    _ScalarMaterials[4]=F;// free energy
    _ScalarMaterials[5]=dFdc;// mu=df/dc->chemical potential
    _ScalarMaterials[6]=d2Fdc2;//dmu/dc
    _ScalarMaterials[7]=Source;
    _ScalarMaterials[8]=dSourcedC;
    _ScalarMaterials[9]=dSourcedMu;
    
}