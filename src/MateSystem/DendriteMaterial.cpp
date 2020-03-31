//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
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

    if(InputParams.size()<7){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for 2d dendrite mate, 7 parameters are required      !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        L,eps,delta,J,theta0,Conduct and eta are expected    !!!   ***\n");
        Msg_AsFem_Exit();
    }

    //*********************************************************************
    //*** Input parameters from your input file should be:
    //*** L      ---> mobility coefficient for Allen-Cahn equation
    //*** eps    ---> interface thickness parameter
    //*** delta  ---> strength of anisotropy
    //*** j      ---> number of angles
    //*** theta0 ---> reference angle
    //*** conduct---> thermal conductivity
    //*** eta    ---> latent heat
    double L=InputParams[0];
    double eps=InputParams[1];
    double delta=InputParams[2];
    int J=int(InputParams[3]);
    double theta0=InputParams[4];
    double Conduct=InputParams[5];
    double eta=InputParams[6];

    // cout<<"L="<<L<<", eps="<<eps<<", delta="<<delta<<", J="<<J<<", theta="<<theta0<<endl;

    double F,dFdphi,d2Fdphi2,d2FdphidT;
    double m,dmdT;

    double phi,T;
    double gradphix,gradphiy;
    const double PI=3.14159265359;
    double K,dKdtheta,d2Kdtheta;// interfacial term and its derivative

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
    const double tol=1.0e-9;
    const double threshold=1.0-tol;
    Vector3d dkdgradphi,ddkdgradphi,v;
    Vector3d dndgradphi;

    n=0.0;
    nsq=gradphix*gradphix+gradphiy*gradphiy;
    if(nsq>tol){
        n=gradphix/sqrt(nsq);
    }

    if(n>threshold){
        n=threshold;
    }
    
    if(n<-threshold){
        n=-threshold;
    }

    theta=acos(n)*sign(gradphiy);

    dthetadn=-sign(gradphiy)/sqrt(1.0-n*n);

    dndgradphi=0.0;
    if(nsq>tol){
        dndgradphi(1)= gradphiy*gradphiy/(nsq*sqrt(nsq));
        dndgradphi(2)=-gradphix*gradphiy/(nsq*sqrt(nsq));
        dndgradphi(3)=0.0;
    }    

    K=eps*(1.0+delta*cos(J*(theta-theta0*PI/180.0)));
    dKdtheta=-eps*delta*J*sin(J*(theta-theta0*PI/180.0));
    d2Kdtheta=-eps*delta*J*J*cos(J*(theta-theta0*PI/180.0));

    dkdgradphi=dKdtheta*dthetadn*dndgradphi;
    ddkdgradphi=d2Kdtheta*dthetadn*dndgradphi;


    // cout<<"dndgradphi="<<dndgradphi(1)<<", "<<dndgradphi(2)<<endl;

    v(1)=-gradphiy;
    v(2)= gradphix;
    v(3)= 0.0;

    

    _ScalarMaterials[0]=L;
    _ScalarMaterials[1]=Conduct;
    _ScalarMaterials[2]=eta;
    _ScalarMaterials[3]=F;
    _ScalarMaterials[4]=dFdphi;
    _ScalarMaterials[5]=d2Fdphi2;
    _ScalarMaterials[6]=d2FdphidT;

    _ScalarMaterials[7]=K;
    _ScalarMaterials[8]=dKdtheta;
    _ScalarMaterials[9]=d2Kdtheta;

    _VectorMaterials[0]=dkdgradphi;
    _VectorMaterials[1]=ddkdgradphi;
    _VectorMaterials[2]=v;
    
}