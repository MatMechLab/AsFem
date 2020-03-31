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

   
    

    double c=gpU[0];// if you use cahn-hilliard,
                    // must make sure the first dof is concentration, not mu!!!
    double tol=1.0e-4;
    if(c<tol) c=tol;
    if(c>1.0-tol) c=1.0-tol;
    _ScalarMaterials[0]=InputParams[0]*c*(1-c);// M=D*c*(1-c)
    _ScalarMaterials[1]=InputParams[0]*(1-2*c);// dM/dc
    _ScalarMaterials[2]=InputParams[1];// Chi
    _ScalarMaterials[3]=InputParams[2];// Kappa
    _ScalarMaterials[4]=c*log(c)+(1-c)*log(1-c)+InputParams[1]*c*(1-c);// free energy
    _ScalarMaterials[5]=log(c)-log(1-c)+InputParams[1]*(1-2*c);// mu=df/dc->chemical potential
    _ScalarMaterials[6]=1.0/c+1.0/(1.0-c)-InputParams[1]*2.0;//dmu/dc
    
}