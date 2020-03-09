//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::DendriteMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for 2d dendrite mate, 5 parameters are required      !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        L, eps, delta, j and K are expected                  !!!   ***\n");
        Msg_AsFem_Exit();
    }

    //*********************************************************************
    //*** Input parameters from your input file should be:
    //*** L    ---> mobility coefficient for Allen-Cahn equation
    //*** eps  ---> interface thickness parameter
    //*** delta---> strength of anisotropy
    //*** j    ---> number of angles
    //*** k    ---> thermal conductivity
    double L=InputParams[0];
    double eps=InputParams[1];
    double delta=InputParams[2];
    int J=int(InputParams[3]);
    double Conduct=InputParams[4];

    double F,dFdphi,d2Fdphi2,d2FdphidT;
    double m,dmdT;

    double phi,T;
    double gradphix,gradphiy;
    const double PI=3.14159265359;
    double K,dK,d2K;// interfacial term and its derivative

    phi=gpU[0];
    gradphix=gpGradU[0](1);
    gradphiy=gpGradU[0](2);

    T=gpU[1];

    m=0.9*atan(10.0*(1-T))/PI;
    dmdT=-9.0/(PI*(1.0+100*(1-T)*(1-T)));

    F=0.25*phi*phi*phi*phi-(0.5-m/3.0)*phi*phi*phi+(0.25-0.5*m)*phi*phi;
    dFdphi=phi*phi*phi-(1.5-m)*phi*phi+(0.5-m)*phi;
    d2Fdphi2=3.0*phi*phi-(3-2*m)*phi+0.5-m;
    d2FdphidT=dmdT*phi*phi-dmdT*phi;

    //******************************************************
    //*** now we do the normal vector's derivative
    //******************************************************
    double theta,dthetadn;
    double n,nsq;
    double tol=1.0e-9;
    Vector3d dkdgradphi,ddkdgradphi,v;
    Vector3d dndgradphi;

    n=0.0;
    nsq=gpGradU[0].normsq();
    if(nsq>tol){
        n=gpGradU[0](1)/gpGradU[0].norm();
    }

    if(n>1.0-tol){
        n=1.0-tol;
    }
    
    if(n<-(1.0-tol)){
        n=-(1.0-tol);
    }

    theta=acos(n)*sign(gpGradU[0](2));

    dthetadn=sign(gpGradU[0](2))/sqrt(1.0-n*n);

    dndgradphi=0.0;
    if(nsq>tol){
        dndgradphi(1)= gradphiy*gradphiy/(gpGradU[0].norm()*gpGradU[0].normsq());
        dndgradphi(2)=-gradphix*gradphiy/(gpGradU[0].norm()*gpGradU[0].normsq());
    }    

    K=eps*(1.0+delta*cos(J*(theta-theta0*PI/180.0)));
    dK=-eps*delta*J*sin(J*(theta-theta0*PI/180.0));
    d2K=-eps*delta*J*J*cos(J*(theta-theta0*PI/180.0));

    dkdgradphi=dK*dthetadn*dndgradphi;
    ddkdgradphi=d2K*dthetadn*dndgradphi;

    v(1)=-gradphiy;
    v(2)= gradphix;
    v(3)= 0.0;

    _ScalarMaterials[0]=L;
    _ScalarMaterials[1]=dFdphi;
    _ScalarMaterials[2]=d2Fdphi2;
    _ScalarMaterials[3]=d2FdphidT;

    _ScalarMaterials[4]=K;
    _ScalarMaterials[5]=dK;
    _ScalarMaterials[6]=d2K;

    _VectorMaterials[0]=dkdgradphi;
    _VectorMaterials[1]=ddkdgradphi;
    _VectorMaterials[2]=v;
    
}