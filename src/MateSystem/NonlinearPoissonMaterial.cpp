//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::NonlinearPoissonMaterial(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<2){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for const poisson mate, two parameters are required  !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        sigma*div(grad(phi))=F, so sigma and F are required  !!!   ***\n");
        Msg_AsFem_Exit();
    }

    //************************
    //*** here the poisson equation is:
    //*** div(sigma*grad(phi))=F
    //**** MateVals[0]-->store sigma
    //**** MateVals[1]-->store dsigma/dphi(for constant case, it is zero)
    //**** MateVals[2]-->store F
    //**** MateVals[3]-->store dF/dphi (for constant case, it is zero)
    if(nDim==1){
        _ScalarMaterials[0]=InputParams[0]*(5.0+sin(gpCoord(1)));// sigma
        _ScalarMaterials[1]=0.0;// dsigma/dphi
        _ScalarMaterials[2]=InputParams[1]*gpU[0];// F
        _ScalarMaterials[3]=InputParams[1];// dF/dphi
    }
    else if(nDim==2){
        // cout<<"work"<<endl;
        _ScalarMaterials[0]=InputParams[0]*(1.25+sin(gpU[0]));// sigma
        _ScalarMaterials[1]=InputParams[0]*cos(gpU[0]);// dsigma/dphi
        _ScalarMaterials[2]=InputParams[1]*(0.0+sin(gpCoord(1)*gpCoord(2)));// F
        _ScalarMaterials[3]=0.0;// dF/dphi
    }
    else if(nDim==3){
        _ScalarMaterials[0]=InputParams[0]*(2.5+sin(gpU[0]));// sigma
        _ScalarMaterials[1]=InputParams[0]*cos(gpU[0]);// dsigma/dphi
        _ScalarMaterials[2]=InputParams[1]*(0.25+sin(gpCoord(1)*gpCoord(2)*gpCoord(3)));// F
        _ScalarMaterials[3]=0.0;// dF/dphi
    }
}