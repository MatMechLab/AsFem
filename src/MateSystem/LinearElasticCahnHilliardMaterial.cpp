//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::LinearElasticCahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
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
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: for linear elastic CH mate, 6 parameters are required !!!  ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        E,nu,D,chi,kappa,cref,omega are expected              !!!  ***\n");
        Msg_AsFem_Exit();
    }

    //*******************************************
    //*** For this model, the DoFs must be: ux,uy,c or ux,uy,uz,c
    //*** Material parameters: E     --->Young's modulus
    //***                      nu    --->Poisson ratio
    //***                      D     --->Diffusion coefficient
    //***                      Chi   --->For two-phase interaction
    //***                      Kappa --->Interfacial energy coefficient
    //***                      Cref  --->Reference concentration
    //***                      Omega --->Partial moler volume(the volume change for 1 mole ions)

    

    const double E    = InputParams[0];
    const double nu   = InputParams[1];
    const double D    = InputParams[2];
    const double Chi  = InputParams[3];
    const double Kappa= InputParams[4];
    const double Cref = InputParams[5];
    const double Omega= InputParams[6];

    _ScalarMaterials[0]=Kappa;

    double F,dFdc,d2Fdc2;
    double M,dMdc;

    // for our elasticity tensor Cijkl
    _Rank4Materials[0].SetToZeros();
    _Rank4Materials[0].SetFromEandNu(E,nu);

    
    RankTwoTensor GradU(0.0),TotalStrain(0.0),ChemicalStrain(0.0),I(0.0);
    I.SetToIdentity();
    double conc;
    if(nDim==2){
        GradU.SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        GradU.SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    // For our total strain
    TotalStrain=(GradU+GradU.Transpose())*0.5;
    // For our chemical strain or eigen-strain
    conc=gpU[nDim];
    const double tol=1.0e-3;
    if(conc<tol    ) conc=tol;
    if(conc>1.0-tol) conc=1.0-tol;
    ChemicalStrain=I*(conc-Cref)*Omega/3.0;
    // For our mechanical strain or the final strain
    _Rank2Materials[0]=TotalStrain-ChemicalStrain;

    // For the mobility term and its derivative
    M=D*conc*(1-conc);
    dMdc=D*(1-2*conc);
    _ScalarMaterials[1]=M;
    _ScalarMaterials[2]=dMdc;



    // this is our stress
    _Rank2Materials[1]=_Rank4Materials[0].DoubleDot(_Rank2Materials[0]);
    // for dStress/dc
    _Rank2Materials[2]=_Rank4Materials[0].DoubleDot(I)*(-Omega/3.0);
    // for our hydrostatic stress SigmaH
    _ScalarMaterials[3]=_Rank2Materials[1].Trace()/3.0;
    // for the derivative of hydrostatic stress w.r.t. c
    _ScalarMaterials[4]=_Rank2Materials[2].Trace()/3.0;

    // for our free energy
    F=conc*log(conc)+(1-conc)*log(1-conc)+Chi*conc*(1-conc)+0.5*_Rank2Materials[1].DoubleDot(_Rank2Materials[0]);
    // for our chemical potential
    dFdc=log(conc)-log(1-conc)+Chi*(1-2*conc)+(-Omega*_ScalarMaterials[3]);
    // for the derivative of chemical potential
    d2Fdc2=1.0/conc+1.0/(1.0-conc)-2*Chi+(-Omega*_ScalarMaterials[4]);

    
    _ScalarMaterials[5]=F;
    _ScalarMaterials[6]=dFdc;
    _ScalarMaterials[7]=d2Fdc2;


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
    _ScalarMaterials[8]=sqrt(1.5*_Rank2Materials[3].DoubleDot(_Rank2Materials[3]));
    if(nDim==2){
        _ScalarMaterials[ 9]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[10]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[11]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[12]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[13]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[14]=_Rank2Materials[0](1,2);//strxy
    }
    else if(nDim==3){
        _ScalarMaterials[ 9]=_Rank2Materials[1](1,1);//sigxx
        _ScalarMaterials[10]=_Rank2Materials[1](2,2);//sigyy
        _ScalarMaterials[11]=_Rank2Materials[1](3,3);//sigzz
        _ScalarMaterials[12]=_Rank2Materials[1](2,3);//sigyz
        _ScalarMaterials[13]=_Rank2Materials[1](1,3);//sigzx
        _ScalarMaterials[14]=_Rank2Materials[1](1,2);//sigxy

        _ScalarMaterials[15]=_Rank2Materials[0](1,1);//strxx
        _ScalarMaterials[16]=_Rank2Materials[0](2,2);//stryy
        _ScalarMaterials[17]=_Rank2Materials[0](3,3);//strzz
        _ScalarMaterials[18]=_Rank2Materials[0](2,3);//stryz
        _ScalarMaterials[19]=_Rank2Materials[0](1,3);//strzx
        _ScalarMaterials[20]=_Rank2Materials[0](1,2);//strxy
    }

}