//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::AnisoPFFractureMaterial(const int &nDim,const double &t,const double &dt,
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

    if(InputParams.size()<10){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for phasefield fracture, 10 parameters are required  !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        Ex,Ey,nux,nuy,Gc,L,viscosity,theta1,2,3  are expected!!!   ***\n");
        Msg_AsFem_Exit();
    }
    // we use the sixth one to indicate wether we use the history one(stager or not)
    int UseHist=1;
    

    //*********************************************
    //*** IMPORTANT!!!
    //*** in this model, d=0--->for undamaged case
    //***                d=1--->for damaged   case
    //**********************************************

    //*******************************************
    //*** material design
    //*** viscosity=ScalarMaterials[0];
    //*** Gc=ScalarMaterials[1];
    //*** L=ScalarMaterials[2];
    //*** H=Hist[0];
    //*******************************************
    const double E1=InputParams[0]; // 
    const double E2=InputParams[1]; //
    const double nu1=InputParams[2]; // 
    const double nu2=InputParams[3]; //

    //Ex,Ey,nux,nuy,Gc,L,viscosity,theta1,2,3

    _ScalarMaterials[0]=InputParams[6];// viscosity
    _ScalarMaterials[1]=InputParams[4];// Gc
    _ScalarMaterials[2]=InputParams[5];// L

    // the euler angle
    double theta1,theta2,theta3;

    theta1=InputParams[7];
    theta2=InputParams[8];
    theta3=InputParams[9];

    RankTwoTensor RotationTensor(0.0);
    RotationTensor.SetFromEulerAngle(theta1,theta2,theta3);

    RankTwoTensor Stress(0.0),StressPos(0.0),StressNeg(0.0);
    RankTwoTensor Strain(0.0),StrainPos(0.0),StrainNeg(0.0);
    
    if(nDim==2){
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2Materials[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    Strain=(_Rank2Materials[2]+_Rank2Materials[2].Transpose())*0.5;
    // our total strain
    _Rank2Materials[0]=Strain;


    RankTwoTensor eigvec(0.0);
    double eigval[3];

    RankFourTensor ProjPos(0.0);
    RankFourTensor I4Sym(RankFourTensor::InitIdentitySymmetric4);
    RankFourTensor ProjNeg(0.0);

    // now we can split the positive and negative stress
    RankTwoTensor I(0.0);

    const double k=5.0e-4; // to avoid the zero stiffness matrix
    double d;
    double Psi,PsiPos,PsiNeg;

    // We use the stress to do the decomposition
    // in this case, we can apply this model to anisotropic case and compressive failure
    // for more details, one is referred to :
    // "https://dukespace.lib.duke.edu/dspace/handle/10161/18247"
    // Yingjie Liu's thesis:
    //     "A Computational Framework for Fracture Modeling in Coupled Field Problems"
    RankFourTensor ElasticityTensor(0.0);
    vector<double> C9(9,0.0);
    double E3,nu3,K1,K2,K3,Mu1,Mu2,Mu3;
    double C11,C12,C13,C22,C23,C33,C44,C55,C66;
    K1=E1/(3*(1-2*nu1));Mu1=E1/(2*(1+nu1));
    
    K2=E2/(3*(1-2*nu2));Mu2=E2/(2*(1+nu2));

    E3=0.5*(E1+E2);
    nu3=0.5*(nu1+nu2);
    
    K3=E3/(3*(1-2*nu3));Mu3=E3/(2*(1+nu3));
    
    C11=(K1+4*Mu1/3.0);
    C12=((K1-2*Mu1/3.0)+(K2-2*Mu2/3.0))*0.5;
    C13=(K1-2*Mu1/3.0+K3-2*Mu3/3.0)*0.5;

    C22=(K2+4*Mu2/3.0);
    C23=(K2-2*Mu2/3.0+K3-2*Mu3/3.0)*0.5;

    C33=K3+4*Mu3/3.0;
    C44=(Mu2+Mu3)*0.5;  //C2323
    C55=(Mu1+Mu3)*0.5;  //C3131
    C66=(Mu1+Mu2)*0.5;  //C1212
    C9[0]=C11;C9[1]=C12;C9[2]=C13;
    C9[3]=C22;C9[4]=C23;
    C9[5]=C33;
    C9[6]=C44;
    C9[7]=C55;
    C9[8]=C66;
    ElasticityTensor.SetFromSymmetric9(C9);
    ElasticityTensor.Rotate(RotationTensor);
    

    Stress  = ElasticityTensor.DoubleDot(Strain);
    ProjPos = Stress.CalcPostiveProjTensor(eigval, eigvec);

    StressPos = ProjPos.DoubleDot(Stress);
    StressNeg = Stress - StressPos;

    // Now we store the final stress:
    d = gpU[nDim];
    const double tol=1.0e-3;
    if (d<    tol) d=tol;
    if (d>1.0-tol) d=1.0-tol;

    // Now the Psi^{+} and Psi^{-} become extremelly easy
    PsiPos = 0.5*StressPos.DoubleDot(Strain);
    PsiNeg = 0.5*StressNeg.DoubleDot(Strain);

    Psi = (1-d)*(1-d)*PsiPos+PsiNeg;

    _Rank2Materials[1]=StressPos*((1-d)*(1-d)+k)+StressNeg;

    // Now its the dStress/dD term
    _Rank2Materials[2] =StressPos*(-2+2*d);

    // For the final jacobian, we can use
    _Rank4Materials[0] =(I4Sym+((1-d)*(1-d)+k-1)*ProjPos).DoubleDot(ElasticityTensor);

    // calculate H, and update the history variable
    _Rank2Materials[3].SetToZeros();// store the dH/dstrain rank2 tensor
    if(PsiPos>=gpHistOld[0]){
        gpHist[0]=PsiPos;
        _Rank2Materials[3]=StressPos;
    }
    else{
        gpHist[0]=gpHistOld[0];
        _Rank2Materials[3].SetToZeros();
    }
    
    //*************************************
    if(UseHist){
        // If we use the history state variable, then the stager solution is used
        _ScalarMaterials[20]=gpHistOld[0];
        _Rank2Materials[5].SetToZeros();
    }
    else{
        // Or we use the fully coupled case
        _ScalarMaterials[20]=gpHist[0];
        _Rank2Materials[5]=_Rank2Materials[3];
    }


    // use the fourth one to calculate the vonMises stress
    _Rank2Materials[4]=_Rank2Materials[1];
    // stress deviator tensor
    _Rank2Materials[4](1,1)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[4](2,2)-=(1.0/3.0)*_Rank2Materials[1].Trace();
    _Rank2Materials[4](3,3)-=(1.0/3.0)*_Rank2Materials[1].Trace();



    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _ScalarMaterials[3]=sqrt(1.5*_Rank2Materials[4].DoubleDot(_Rank2Materials[4]));
    // hydrostatic stress
    _ScalarMaterials[4]=_Rank2Materials[1].Trace()/3.0;
    if(nDim==2){
        _ScalarMaterials[ 5]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[ 6]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[ 7]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[ 8]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[ 9]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[10]=_Rank2Materials[0](1,2);//strxy

        _ScalarMaterials[11]=Psi;
        _ScalarMaterials[12]=PsiPos;
        _ScalarMaterials[13]=PsiNeg;
    }
    else if(nDim==3){
        _ScalarMaterials[ 5]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[ 6]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[ 7]=_Rank2Materials[1](3,3);//sigzz
        _ScalarMaterials[ 8]=_Rank2Materials[1](2,3);//sigyz
        _ScalarMaterials[ 9]=_Rank2Materials[1](1,3);//sigzx
        _ScalarMaterials[10]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[11]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[12]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[13]=_Rank2Materials[0](3,3);//strzz
        _ScalarMaterials[14]=_Rank2Materials[0](2,3);//stryz
        _ScalarMaterials[15]=_Rank2Materials[0](1,3);//strzx
        _ScalarMaterials[16]=_Rank2Materials[0](1,2);//strxy

        _ScalarMaterials[17]=Psi;
        _ScalarMaterials[18]=PsiPos;
        _ScalarMaterials[19]=PsiNeg;
    }

}