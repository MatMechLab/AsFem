//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::UserMaterial1(const int &nDim,const double &t,const double &dt,
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
    if(InputParams[0]){}

    // if(InputParams.size()<2){
    //     PetscPrintf(PETSC_COMM_WORLD,"*** Error: for user material1, two parameters are required  !!!   ***\n");
    //     PetscPrintf(PETSC_COMM_WORLD,"***        sigma*div(grad(phi))=F, so sigma and F are required  !!!   ***\n");
    //     Msg_AsFem_Exit();
    // }

    //************************
    //*** here the nonlinear poisson equation is:
    //*** -div((1+u^2)*grad(phi))=F
    //**** MateVals[0]-->store sigma
    //**** MateVals[1]-->store dsigma/dphi(for constant case, it is zero)
    //**** MateVals[2]-->store F
    //**** MateVals[3]-->store dF/dphi (for constant case, it is zero)

    _ScalarMaterials[0]=1.0;// sigma
    _ScalarMaterials[1]=0.0;// dsigma/dphi
    _ScalarMaterials[2]=gpCoord(1)*sin(gpCoord(2));// F=x*sin(y)
    _ScalarMaterials[3]=0.0;// dF/dphi

}