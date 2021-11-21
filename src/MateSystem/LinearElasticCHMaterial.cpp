//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.11.21
//+++ Purpose: Implement the calculation of mechanically coupled 
//+++          CahnHilliard equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/LinearElasticCHMaterial.h"


void LinearElasticCHMaterial::InitMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,Materials &Mate){
    // Here we do not consider any initial internal strains, stress
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}

    Mate.Rank2Materials("stress").SetToZeros();
}

//*********************************************************
void LinearElasticCHMaterial::ComputeDeformationGradientTensor(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &F){
    // here we calculate the deformation gradient tensor as well as euler-lagrange strain tensor
    // DoFs:
    //   1) concentration
    //   2) mu
    //   3) ux
    //   4) uy
    //   5) uz
    if(elmtinfo.nDim==1){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==2){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    else if(elmtinfo.nDim==3){
        _GradU.SetFromGradU(elmtsoln.gpGradU[3],elmtsoln.gpGradU[4],elmtsoln.gpGradU[5]);
    }
    F=_GradU;// for small strain case, we set F to be GradU !
}
//****************************************************************
void LinearElasticCHMaterial::CalcFreeEnergyAndDerivatives(const vector<double> &InputParams, const double &c,double &df,double &d2f){
    // f=(x-xa)^2*(x-xb)^2;
    if(InputParams.size()){}// in this case, we do not use the xa, xb from input file
    double height=InputParams[2-1];//energy barrier height
    const double ca=0.1;
    const double cb=0.9;
    df=height*2*(c-ca)*(c-cb)*(2*c-ca-cb);
    d2f=height*2*(ca*ca+4*ca*cb+cb*cb-6*ca*c-6*cb*c+6*c*c);
}
//****************************************************************
void LinearElasticCHMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const LocalElmtSolution &soln,const RankTwoTensor &F,RankTwoTensor &Strain,RankTwoTensor &Stress,RankFourTensor &Jacobian){
    
    double E,nu,C,Omega,C0;
    double dFdc,d2Fdc2;

    // take parameters from InputParams
    E=InputParams[4-1];
    nu=InputParams[5-1];
    Omega=InputParams[6-1];
    C0=InputParams[7-1];

    C=soln.gpU[1];

    // calculate the chemical potential and its derivative
    CalcFreeEnergyAndDerivatives(InputParams,C,dFdc,d2Fdc2);
    
    // for different strains(coupling term) 
    _TotalStrain=(F+F.Transpose())*0.5;// the total strain, here F is GradU
    _I.SetToIdentity();
    // for the concentration induced eigen strain 
    _EigenStrain=_I*(C-C0)*Omega/3.0;
    _dEigenStraindC=_I*(Omega/3.0);
    // for elastic strain
    _MechStrain=_TotalStrain-_EigenStrain;
    _dMechStraindC=_dEigenStraindC*-1.0;
    Strain=_MechStrain;

    // set up the final jacobian, in this case, its simple linear elasticity tensor
    Jacobian.SetToZeros();
    Jacobian.SetFromEandNu(E,nu);

    // the total free energy density is:
    // psi=psic+psie
    // psic=(c-ca)^2(c-cb)^2
    // psie=0.5*(strain-eigenstrain):C:(strain-eigenstrain)
    Stress=Jacobian.DoubleDot(_MechStrain);
    _dStressdC=Jacobian.DoubleDot(_dMechStraindC);

    // now we calculate the final chemical potentials
    // psie=0.5*sigma:strain
    _Mu=dFdc+0.5*_dStressdC.DoubleDot(_MechStrain)+0.5*Stress.DoubleDot(_dMechStraindC);
    _dMudC=d2Fdc2+_dStressdC.DoubleDot(_dMechStraindC);

    _dMudStrain=_dStressdC;
}
//***************************************************************
void LinearElasticCHMaterial::ComputeMaterialProperties(const vector<double> &InputParams,const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,
                                           const Materials &MateOld,Materials &Mate) {

    //*********************************************************
    //*** get rid of unused warnings
    //*********************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    if(InputParams.size()<7){
        MessagePrinter::PrintErrorTxt("for the LinearElasticCHMaterial, seven parameters are required, you need to give: D, barrier height, kappa, E, nu, Omega, c0");
        MessagePrinter::AsFem_Exit();
    }
    Mate.ScalarMaterials("M")=InputParams[1-1];
    Mate.ScalarMaterials("dMdC")=0.0;
    Mate.ScalarMaterials("Kappa")=InputParams[3-1];

    ComputeDeformationGradientTensor(elmtinfo,elmtsoln,_F);
    ComputeConstitutiveLaws(InputParams,elmtsoln,_F,_Strain,_Stress,_Jac);


    _devStress=_Stress.Dev();
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_devStress.DoubleDot(_devStress));
    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jac;

    Mate.ScalarMaterials("dFdC")=_Mu;
    Mate.ScalarMaterials("d2FdC2")=_dMudC;
    Mate.ScalarMaterials("HyStress")=_Stress.Trace()/3.0;
    Mate.Rank2Materials("dMudStrain")=_dMudStrain;
    Mate.Rank2Materials("dStressdC")=_dStressdC;

}
