//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::TensorCahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for tensor CH mate, seven parameters are required    !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        D_xx,D_yy,D_zz[zy,zx,xy],chi,Kappa_xx,yy,zz[zy,zx,xy]!!!   ***\n");
        Msg_AsFem_Exit();
    }

   
    

    double c=gpU[0];// if you use cahn-hilliard,
                    // must make sure the first dof is concentration, not mu!!!
    double tol=1.0e-3;
    if(c<tol) c=tol;
    if(c>1.0-tol) c=1.0-tol;

    double dxx,dyy,dzz,dzy,dzx,dxy;
    double kxx,kyy,kzz,kzy,kzx,kxy;
    double chi;
    dxx=0.0;dyy=0.0;dzz=0.0;dzy=0.0;dzx=0.0;dxy=0.0;
    kxx=0.0;kyy=0.0;kzz=0.0;kzy=0.0;kzx=0.0;kxy=0.0;
    if(InputParams.size()==7){

        dxx=InputParams[0];
        dyy=InputParams[1];
        dzz=InputParams[2];
        dzy=0.0;dzx=0.0;dxy=0.0;
        
        chi=InputParams[3];
        
        _ScalarMaterials[0]=chi;// Chi
        _ScalarMaterials[1]=c*log(c)+(1-c)*log(1-c)+chi*c*(1-c);// free energy
        _ScalarMaterials[2]=log(c)-log(1-c)+chi*(1-2*c);        // mu=df/dc->chemical potential
        _ScalarMaterials[3]=1.0/c+1.0/(1.0-c)-chi*2.0;          //dmu/dc

        kxx=InputParams[4];
        kyy=InputParams[5];
        kzz=InputParams[6];
        kzy=0.0;kzx=0.0;kxy=0.0;
    }
    else if(InputParams.size()>=13){
        dxx=InputParams[ 0];
        dyy=InputParams[ 1];
        dzz=InputParams[ 2];
        dzy=InputParams[ 3];
        dzx=InputParams[ 4];
        dxy=InputParams[ 5];

        chi=InputParams[6];
        _ScalarMaterials[0]=chi;// Chi
        _ScalarMaterials[1]=c*log(c)+(1-c)*log(1-c)+chi*c*(1-c);// free energy
        _ScalarMaterials[2]=log(c)-log(1-c)+chi*(1-2*c);// mu=df/dc->chemical potential
        _ScalarMaterials[3]=1.0/c+1.0/(1.0-c)-chi*2.0;//dmu/dc

        kxx=InputParams[ 7];
        kyy=InputParams[ 8];
        kzz=InputParams[ 9];
        kzy=InputParams[10];
        kzx=InputParams[11];
        kxy=InputParams[12];
    }

    _Rank2Materials[0](1,1)=dxx;_Rank2Materials[0](1,2)=dxy;_Rank2Materials[0](1,3)=dzx;
    _Rank2Materials[0](2,1)=dxy;_Rank2Materials[0](2,2)=dyy;_Rank2Materials[0](2,3)=dzy;
    _Rank2Materials[0](3,1)=dzx;_Rank2Materials[0](3,2)=dzy;_Rank2Materials[0](3,3)=dzz;

    _Rank2Materials[1]=_Rank2Materials[0];
    _Rank2Materials[0]*=c*(1.0-c);
    _Rank2Materials[1]*=(1.0-2.0*c);

    // cout<<_Rank2Materials[0](1,1)<<" "<<_Rank2Materials[0](1,2)<<" "<<_Rank2Materials[0](1,3)<<endl;
    // cout<<_Rank2Materials[0](2,1)<<" "<<_Rank2Materials[0](2,2)<<" "<<_Rank2Materials[0](2,3)<<endl;
    // cout<<_Rank2Materials[0](3,1)<<" "<<_Rank2Materials[0](3,2)<<" "<<_Rank2Materials[0](3,3)<<endl<<endl;

    // cout<<_Rank2Materials[1](1,1)<<" "<<_Rank2Materials[1](1,2)<<" "<<_Rank2Materials[1](1,3)<<endl;
    // cout<<_Rank2Materials[1](2,1)<<" "<<_Rank2Materials[1](2,2)<<" "<<_Rank2Materials[1](2,3)<<endl;
    // cout<<_Rank2Materials[1](3,1)<<" "<<_Rank2Materials[1](3,2)<<" "<<_Rank2Materials[1](3,3)<<endl<<endl;
    // cout<<"****************"<<endl<<endl;

    _Rank2Materials[2](1,1)=kxx;_Rank2Materials[2](1,2)=kxy;_Rank2Materials[2](1,3)=kzx;
    _Rank2Materials[2](2,1)=kxy;_Rank2Materials[2](2,2)=kyy;_Rank2Materials[2](2,3)=kzy;
    _Rank2Materials[2](3,1)=kzx;_Rank2Materials[2](3,2)=kzy;_Rank2Materials[2](3,3)=kzz;
    
}