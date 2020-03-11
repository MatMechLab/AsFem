//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::CurrentThermalMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for const poisson mate, 5 parameters are required    !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        C, lambda, R, I and rho are expected                 !!!   ***\n");
        Msg_AsFem_Exit();
    }

    //************************
    //*** here the poisson equation is:
    //*** C      ---> heat capacity
    //*** lambda ---> thermal diffusivity
    //*** R      ---> resistivity
    //*** I      ---> applied current
    //*** rho    ---> density

    double C,lambda,R,I,rho;

    C=InputParams[0];
    lambda=InputParams[1];
    R=InputParams[2];
    I=InputParams[3];
    rho=InputParams[4];
    
    // for linear case
    double Alpha,dAlphadT,F,dFdT,Q;

    // for the joule heat
    Q=I*I*R*t;

    // for the thermal conductivity
    Alpha=lambda/(rho*C);
    dAlphadT=0.0;

    // for thermal source term
    F=Q/(rho*C);
    dFdT=0.0;


    // store them into our material database
    _ScalarMaterials[0]=Alpha;
    _ScalarMaterials[1]=dAlphadT;
    _ScalarMaterials[2]=F;
    _ScalarMaterials[3]=dFdT;

    
}