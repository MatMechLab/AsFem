//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.11.21
//+++ Purpose: implement the residual and jacobian for general
//+++          mechanically coupled CahnHilliard equation
//+++          1) dc/dt=div(M*grad(Mu))
//+++          2) mu   =dF/dc-kappa*grad(c)
//+++          3) div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/MechanicsCahnHilliardElmt.h"

void MechanicsCahnHilliardElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) {
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(elmtinfo,soln,shp,Mate,MateOld,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(elmtinfo,ctan,soln,shp,Mate,MateOld,localK);
    }
    else if(calctype==FECalcType::Projection){
        ComputeProjection(elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in MechanicsCahnHilliardElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*****************************************************************************
void MechanicsCahnHilliardElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    //******************************************************************
    // calculate the residual contribution of CahnHilliard equation
    //******************************************************************
    // For R_c
    M=Mate.ScalarMaterials("M");
    dFdC=Mate.ScalarMaterials("dFdC");
    kappa=Mate.ScalarMaterials("Kappa");
    localR(1)=soln.gpV[1]*shp.test
        +M*soln.gpGradU[2]*shp.grad_test;
    // For R_mu
    localR(2)=soln.gpU[2]*shp.test-dFdC*shp.test-kappa*soln.gpGradU[1]*shp.grad_test;
    //***************************************************
    // For mechanics part
    Stress=Mate.Rank2Materials("stress")-MateOld.Rank2Materials("stress");
    localR(3)=Stress.IthRow(1)*shp.grad_test;
    if(elmtinfo.nDim>=2){
        localR(4)=Stress.IthRow(2)*shp.grad_test;
        if(elmtinfo.nDim==3){
            localR(5)=Stress.IthRow(3)*shp.grad_test;
        }
    }
}
//******************************************************************************
void MechanicsCahnHilliardElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    M=Mate.ScalarMaterials("M");
    dMdC=Mate.ScalarMaterials("dMdC");
    d2FdC2=Mate.ScalarMaterials("d2FdC2");
    kappa=Mate.ScalarMaterials("Kappa");
    dStressdC=Mate.Rank2Materials("dStressdC");
    dMudStrain=Mate.Rank2Materials("dMudStrain");

    //**********************************************
    //*** contribution of C
    //**********************************************
    // K_c,c
    localK(1,1)=shp.trial*shp.test*ctan[1]
        +dMdC*shp.trial*soln.gpGradU[2]*shp.grad_test*ctan[0];
    // K_c,mu
    localK(1,2)=M*shp.grad_trial*shp.grad_test*ctan[0];
    // K_c,ui
    localK(1,3)=0.0;                      // k_c,ux
    if(elmtinfo.nDim>=2) localK(1,4)=0.0; // k_c,uy
    if(elmtinfo.nDim==3) localK(1,5)=0.0; // k_c,uz
    //*********************************************
    //*** contribution of mu
    //*********************************************
    // K_mu,c
    localK(2,1)=-d2FdC2*shp.trial*shp.test*ctan[0]
        -kappa*shp.grad_trial*shp.grad_test*ctan[0];
    // K_mu,mu
    localK(2,2)=shp.trial*shp.test*ctan[0];
    // for the mechanical coupling
    valx=0.0;valy=0.0;valz=0.0;
    for(int i=1;i<=3;i++){
        valx+=0.5*(dMudStrain(1,i)+dMudStrain(i,1))*shp.grad_trial(i);
        valy+=0.5*(dMudStrain(2,i)+dMudStrain(i,2))*shp.grad_trial(i);
        valz+=0.5*(dMudStrain(3,i)+dMudStrain(i,3))*shp.grad_trial(i);
    }
    // K_mu,ux
    localK(2,3)=-valx*shp.test*ctan[0];                     // K_mu,ux
    if(elmtinfo.nDim>=2) localK(2,4)=-valy*shp.test*ctan[0];// K_mu,uy
    if(elmtinfo.nDim==3) localK(2,5)=-valz*shp.test*ctan[0];// K_mu,uz
    
    //*********************************************
    //*** for stress equilibrium equation
    //*********************************************
    // K_ux,ux
    localK(3,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,1,shp.grad_test,shp.grad_trial)*ctan[0];
    // K_ux,c
    localK(3,1)=dStressdC.IthRow(1)*shp.grad_test*shp.trial*ctan[0];
    // K_ux,mu
    localK(3,2)=0.0;
    if(elmtinfo.nDim>=2){
        // K_ux,uy
        localK(3,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,2,shp.grad_test,shp.grad_trial)*ctan[0];
        //*************************
        // K_uy,c
        localK(4,1)=dStressdC.IthRow(2)*shp.grad_test*shp.trial*ctan[0];
        // K_uy,mu
        localK(4,2)=0.0;
        // K_uy,ux
        localK(4,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,1,shp.grad_test,shp.grad_trial)*ctan[0];
        // K_uy,uy
        localK(4,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,2,shp.grad_test,shp.grad_trial)*ctan[0];
        if(elmtinfo.nDim==3){
            // K_ux,uz
            localK(3,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uy,uz
            localK(4,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,c
            localK(5,1)=dStressdC.IthRow(3)*shp.grad_test*shp.trial*ctan[0];
            // K_uz,mu
            localK(5,2)=0.0;
            // K_uz,ux
            localK(5,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,1,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uy
            localK(5,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,2,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uz
            localK(5,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,3,shp.grad_test,shp.grad_trial)*ctan[0];
        }
    }
}
//*************************************************
void MechanicsCahnHilliardElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}

}
