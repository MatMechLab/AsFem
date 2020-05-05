//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::SaintVenantMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for saint-venant mate, two parameters are required   !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        E and nu is expected for elastic problem             !!!   ***\n");
        Msg_AsFem_Exit();
    }


    const double EE=InputParams[0];
    const double nu=InputParams[1];
    //***********************************************
    //*** the related formulas are taken from:
    //*** https://en.wikipedia.org/wiki/Hyperelastic_material
    //************************************************
    //**** For our Cijkl elasticity tensor
    _Rank4Materials[0].SetToZeros();
    _Rank4Materials[0].SetFromEandNu(EE,nu);


    //*******************************************
    //*** now we calculate the finite strain
    //*******************************************
    _Rank2Materials[2].SetToZeros();
    RankTwoTensor F(0.0),I(0.0),C(0.0),E(0.0);
    
    I.SetToIdentity();
    
    if(nDim==2){
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }

   
    F=_Rank2Materials[2]+I;// F=I+U_{i,j}
    C=F.Transpose()*F;//C=F^tF

    E=(C-I)*0.5;


    //******************************
    //*** our finite strain
    //******************************
    _Rank2Materials[0]=E;

    //******************************
    //*** our stress
    //******************************
    _Rank2Materials[1]=_Rank4Materials[0].DoubleDot(E);


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
    }

}